# -*- coding: utf-8 -*-
## 2018/10/26: copy and modified classification script from hpc

import sys
import os
from sklearn.model_selection import StratifiedKFold, KFold
import pandas as pd
import numpy as np
import time
from functools import reduce

import sklearn
import xgboost
#print("sklearn version: ", sklearn.__version__)
#print("xgboost version: ", xgboost.__version__)

## import custom functions
codepath = "~/EVE/eve/models/"
sys.path.append(codepath)
from XGBoost_utils import unit_train, get_RFE_sizes

## remove numpy - sklearn bug warning ##
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
from sklearn.exceptions import DataConversionWarning
warnings.filterwarnings(action='ignore', category=DataConversionWarning)
## https://stackoverflow.com/questions/49545947/sklearn-deprecationwarning-truth-value-of-an-array

start_time = time.time()

#################
## User Inputs ##
#################
runSpec = {}
with open(sys.argv[2], "r") as f:
  for line in f:
    if "$" in line:
      title = line.replace("$","").strip()

    if "[1]" in line:
      val = line.replace("[1] ","").strip()
      val = val.replace('"', '')

      runSpec[title] = val

seed = int(sys.argv[1])
cv_id_curr = int(sys.argv[3])

label = runSpec["label_name"]
featureList = [] ## input list of features. Use all features if it is left empty.
n_cv = int(runSpec["num_CV"]) ## number of CV

RFE_step = runSpec["RFE_step"].split()
RFE_step = [float(num) for num in RFE_step]
RFE_criteria = runSpec["RFE_criteria"] ## string, gain or freq to use in RFE

if runSpec["surv_col"] == "NA":
  col_survival = []
else:
  col_survival = runSpec["surv_col"] # "PFS" # can be []

if runSpec["event_col"] == "NA":
  col_event = []
else:
  col_event = runSpec["event_col"] #"PFS.event" # can be []

## =============================================================================
## Check whether files already exist, if so, stop the process.
## =============================================================================
outP = runSpec["outP"] ## get absolute path
outP = os.path.expanduser(outP)

out_vimp = outP +"/xgb_vimp_"+label+"_seed"+str(seed)+"_cv"+str(cv_id_curr)+".csv"
out_preval = outP +"/xgb_prevalidation_"+label+"_seed"+str(seed)+"_cv"+str(cv_id_curr)+".csv"
out_grid = outP +"/xgb_gridsearch_"+label+"_seed"+str(seed)+"_cv"+str(cv_id_curr)+".csv"

if os.path.exists(out_vimp):
  #print("File already exist, quitting")
  sys.exit("File already exist, quitting ...")

## =============================================================================
## advanced setting 
## =============================================================================


## read data
df = pd.read_csv(runSpec["project_home"]+ "/" + runSpec["training_data"])

## collect sample ID info for pre-validation
if runSpec["sample_ID"] == "NA":
    sampleID = "RowIndex"
    df_sampleid = df.index.values
else:
    sampleID = runSpec["sample_ID"]
    df_sampleid = df[sampleID]

## preproc x & y
## whether to calculate hazard ratio
if col_survival and col_event:
    ## if user provided both survival time and event flags
    y = df.loc[:,[col_survival, col_event]]
    y.columns = ["col_surv", "col_event"]
    HR_calc = True
else:
    ## if user only provide label 
    y = df.loc[:,label]
    y.columns = ["lb"]
    HR_calc = False

if not featureList:
    print("Use all features")
    col2drop = [label, sampleID, "USUBJID", "OS", "OS.event", "PFS", "PFS.event", runSpec["weight_col"]]
    idx = df.columns.isin(col2drop)
    x = df.iloc[:, ~idx] ## remove those columns that are not features
else:
    idx = df.columns.isin(featureList)
    x = df.iloc[:,idx] ## select custom selected features

## handel RFE
sizes = get_RFE_sizes(RFE_step, x.shape[1])
print("RFE steps:", sizes)

## objective and evaluation functions
if label in ["OS", "PFS"]:
    obj = "survival:cox"
    evalm = "cox-nloglik"
    
    ## re-modify event flag for XGBoost
    event_flag = np.array([i if i == 1 else -1 for i in df.loc[(label+".event")]])
    y = y * event_flag

elif df[label].nunique() > 2: ## if more than 2 classes  
    obj = "multi:softmax"
    evalm = "mlogloss" ## merror

elif df[label].nunique() == 2: ## if binary classes
    obj = "binary:logistic" 
    evalm = "logloss" ## "auc""
    
else:
    print("Label contains only 1 class.")
    
## xgboost hyper-params
xgb_params = {
    "n_estimators": 3000,
    "learning_rate": 0.01, # 0.3, [0,1], usually 0.01 - 0.2
    "gamma": 0, # 0, [0,∞], minimum loss reduction required to make a further partition on a leaf node of the tree. The larger, the more conservative the algorithm will be.
    "max_depth": 4, # 6, 0 means no limit. [0,∞]
    "min_child_weight": 1, # 1, [0,∞], If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. In linear regression mode, this simply corresponds to minimum number of instances needed to be in each node. The larger, the more conservative the algorithm will be.
    "max_delta_step": 0, # 0, [0,∞], Usually this parameter is not needed, but it might help in logistic regression when class is extremely imbalanced. Set it to value of 1-10 might help control the update
    "subsample": 1, # 1, [0,1], subsample ratio of the training instance. Setting it to 0.5 means that XGBoost randomly collected half of the data instances to grow trees and this will prevent overfitting.
    "colsample_bytree": 0.8, # 1, [0,1], subsample ratio of columns when constructing each tree
    "reg_lambda": 1, # 1, L2 regularization term on weights, increase this value will make model more conservative
    "reg_alpha": 0, # 0, L1 regularization term on weights, increase this value will make model more conservative
    #"scale_pos_weight": float(runSpec["scale_pos_weight"]), #581/219, #1,  sum(negative cases) / sum(positive cases)
    "silent": 1,
    "objective": obj,
    "njobs": -1,
    "random_state": seed
}

####################
## training start ##
####################

## check if LOOCV
if df.shape[0] == n_cv:
  ss = KFold(n_splits=n_cv, shuffle=True, random_state=seed*2)
else:
  ss = StratifiedKFold(n_splits=n_cv, shuffle=True, random_state=seed*2)

## set up initial xgboost
if label=="OS" or label=="PFS":
  xgbc = xgboost.sklearn.XGBRegressor(**xgb_params)
  ss.get_n_splits(x)
else:
  xgbc = xgboost.sklearn.XGBClassifier(**xgb_params)
  ss.get_n_splits(x, y)
    
## initialize empty output dataframes
df_vimp = pd.DataFrame()
#df_error = pd.DataFrame()
df_grid = pd.DataFrame()
df_prevalid = pd.DataFrame()

cv_id = 1
for train_index, test_index in ss.split(x, df[label]):
  #print("TRAIN:", train_index, "TEST:", test_index)
  
  ## get weights
  if runSpec["weight_col"] == "NA":
    weights = None
  else:
    weights = df.loc[train_index, runSpec["weight_col"]]
  
  if (cv_id_curr == 0):
    print("Fold number:", cv_id)
    X_train, X_test = x.iloc[train_index], x.iloc[test_index]
    Y_train, Y_test = y.iloc[train_index], y.iloc[test_index]
    
    df_vimp_tmp, df_grid_tmp, df_prevalid_tmp \
        = unit_train(xgbc, X_train, Y_train, X_test, Y_test, 
                     sizes, evalm, HR_calc, RFE_criteria, weights)
    
    df_vimp_tmp["cv"] = cv_id
    #df_error_tmp["cv"] = cv_id
    df_grid_tmp["cv"] = cv_id
    df_prevalid_tmp["cv"] = cv_id
    
    df_vimp = df_vimp.append(df_vimp_tmp)
    #df_error = df_error.append(df_error_tmp)
    df_grid = df_grid.append(df_grid_tmp)
    df_prevalid = df_prevalid.append(df_prevalid_tmp)
  
  elif (cv_id_curr == cv_id):
    print("Fold number:", cv_id)
    X_train, X_test = x.iloc[train_index], x.iloc[test_index]
    Y_train, Y_test = y.iloc[train_index], y.iloc[test_index]
    
    df_vimp_tmp, df_grid_tmp, df_prevalid_tmp \
        = unit_train(xgbc, X_train, Y_train, X_test, Y_test, 
                     sizes, evalm, HR_calc, RFE_criteria, weights)
  
    df_vimp_tmp["cv"] = cv_id
    #df_error_tmp["cv"] = cv_id
    df_grid_tmp["cv"] = cv_id
    df_prevalid_tmp["cv"] = cv_id
    
    df_vimp = df_vimp.append(df_vimp_tmp)
    #df_error = df_error.append(df_error_tmp)
    df_grid = df_grid.append(df_grid_tmp)
    df_prevalid = df_prevalid.append(df_prevalid_tmp)
  
  cv_id = cv_id + 1

## measure duration of the unit code run
end_time = time.time()
time_taken = end_time - start_time
print("time taken:", time_taken)

## process the pre-validation
# 1. merge sample id to the prevalidation data based on index
df_prevalid2 = pd.merge(df_prevalid, pd.DataFrame({sampleID: df_sampleid}), 
    left_index=True, right_index=True)

## save output to csv
if not os.path.exists(outP):
  print("create results path:", outP)
  os.makedirs(outP)
    
#out_vimp = outP +"/xgb_vimp_"+label+"_seed"+str(seed)+"_cv"+str(cv_id_curr)+".csv"
#out_preval = outP +"/xgb_prevalidation_"+label+"_seed"+str(seed)+"_cv"+str(cv_id_curr)+".csv"
#out_grid = outP +"/xgb_gridsearch_"+label+"_seed"+str(seed)+"_cv"+str(cv_id_curr)+".csv"

df_vimp.to_csv(out_vimp, index=True)
df_prevalid2.to_csv(out_preval, index=False)
df_grid.to_csv(out_grid, index=False)
