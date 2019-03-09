## 2018/10/03: implement pred_proba and calibration
## 2018/09/20: implement pre-validation, use 'Gain' to do RFE (added an option for user to select what to use in RFE)

import sys
from sklearn.model_selection import GridSearchCV, KFold, RandomizedSearchCV
from sklearn.externals import joblib 
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score, roc_auc_score
from sklearn.calibration import CalibratedClassifierCV
#from lifelines import CoxPHFitter
import pandas as pd
import numpy as np

## GridSearch for max_depth and min_child_weight 
def grid_search(xgbc, X_train, Y_train, X_test, Y_test, evalm, weights=None):
    """ perform Randomized Search
    n_iter trades of between performance and speed.
    n_iter large will be just like gridsearch, but will have the best performance.
    """
    print("grid search starts ...")

    parameters = {'max_depth':[1,2,3,4,5], 
                  #'learning_rate':[0.1, 0.05],
                  #'min_child_weight':[1,3,5],
                  #"n_estimators":[100, 500, 1000, 1500, 2000],
                  "max_delta_step":[0,1,3,5],
                  "reg_alpha":[0,1,3,5]
                 }

    fit_params={"early_stopping_rounds":100, 
                "sample_weight": weights,
                "eval_metric" : evalm, 
                "eval_set" : [(X_test, Y_test)],
                "verbose": False}
    
    cv_sets = KFold(n_splits=3, shuffle=True, random_state=42)
    
    if evalm == "cox-nloglik":
        parameters["n_estimators"] = [100, 300, 500]
        grid = RandomizedSearchCV(xgbc, parameters, cv=cv_sets, 
                        fit_params = fit_params, 
                        n_jobs=-1, error_score = 0, n_iter = 20)
    elif evalm == "rmse":
        grid = RandomizedSearchCV(xgbc, parameters, cv=cv_sets, 
                        fit_params = fit_params, 
                        n_jobs=-1, scoring = "neg_mean_squared_error",
                        error_score = 0, n_iter = 20)
    else:
        grid = RandomizedSearchCV(xgbc, parameters, cv=cv_sets, 
                            fit_params = fit_params,
                            scoring = 'f1_macro', n_jobs=-1,
                            error_score = 0, n_iter = 20)

    with joblib.parallel_backend('threading'):
        grid.fit(X_train, Y_train)
        print(joblib.effective_n_jobs(), "effective CPU used")
    
    print("The best parameters are %s with a score of %0.2f" % (grid.best_params_, grid.best_score_))
    return(grid)

# def get_hr(Y_test_surv, pred):
#     """Calculates hazard ratio between 2 classes only
#     
#     Parameters
#     ----------
#     Y_test_surv: dataframe 
#         It has to contains two columns, 
#         one is survival time (PFS/OS), one is event flag (1/0)
#         colname has to be col_surv and col_event respectively
#     pred: np.array or pd.series
#         This is the prediciton results (array) for each patients
#     
#     Returns
#     ---------
#     hr : float
#         hazard ratio for the class 1 vs class 0
#         
#     """
#     df = Y_test_surv.copy()
#     df['pred'] = pred
#     cox = CoxPHFitter()
#     cox.fit(df, 'col_surv', 'col_event')
#     hr = np.exp(cox.hazards_)['pred'].values
#     
#     return(np.float(hr))

def metrics(y_true, y_pred, y_pred_prob = None):
    '''Calculate some model evaluation metrics for pre-validation
    
    Parameters
    ----------
    y_true: pd.Series or np.array
        true labels
    y_pred: pd.Series or np.array
        predicted labels
    y_pred_prob: pd.Series or np.array
        probabilities of the predicted labels
    
    Return
    ----------
    df_dict: dictionary
        columns might include: 
            f1 score for each label,
            accuracy, AUC, ...
        
    '''

    lbs = y_true.unique()
    lbs.sort()
    #recall = recall_score(y_true, y_pred, labels=lbs, average=None) ## specificity
    #precision = precision_score(y_true, y_pred, labels=lbs, average=None)
    acc = accuracy_score(y_true, y_pred)
    
    f1 = f1_score(y_true, y_pred, labels=lbs, average=None)
    df_dict = {}
    for i, lb in enumerate(lbs):
        df_dict["f1_"+str(lb)] = f1[i]
    df_dict["acc"] = acc
    
    if y_pred_prob:
      auc = roc_auc_score(y_true, y_pred_prob)
      df_dict["rocauc"] = auc
    
    #df_dict = pd.Series(df_dict).to_frame().T
    
    return df_dict

def confM(xgbc, X_test, Y_test, evalm, xgbc_cali=None, Y_test_surv=None):
    ## ToDo: this might cause problem, if Y_test does not contain all labels
    lbs = Y_test.unique() 
    lbs.sort()
    Y_pred = xgbc.predict(X_test)
    
    #df_dict = {}
    #df_dict[evalm] = xgbc.evals_result_['validation_0'][evalm][-1]
    
    if evalm in ["auc", "merror", "logloss", "mlogloss"]:
      Y_pred_cali = xgbc_cali.predict(X_test)
      Y_pred_prob = xgbc.predict_proba(X_test) ## might not work for regression
      Y_pred_prob_cali = xgbc_cali.predict_proba(X_test) ## calibrated outcomes
      
      ## for pre-validation
      df_prevalid = Y_test.to_frame().copy()
      df_prevalid['pred'] = Y_pred
      df_prevalid['pred_cali'] = Y_pred_cali
      
      cm = confusion_matrix(Y_test, Y_pred, labels=lbs)
      
      for i, lb in enumerate(lbs):
          #tot_count = sum(cm[i,])
          #cm[i,i] = 0 ## make the diaganol = 0
          #error_count = sum(cm[i,])
          #df_dict[str(lb)] = int(error_count)
          #df_dict["n_"+str(lb)] = int(tot_count)
          df_prevalid['predprob_'+str(lb)] = Y_pred_prob[:, i]
          df_prevalid['predprob_cali_'+str(lb)] = Y_pred_prob_cali[:, i]
      
      #df_dict_tmp = metrics(Y_test, Y_pred) # Y_pred_prob
      #df_dict.update(df_dict_tmp)
    
    elif evalm == "cox-nloglik":
        df_prevalid = Y_test_surv.copy()
        df_prevalid['pred'] = Y_pred
    
    elif evalm == "rmse":
        df_prevalid = Y_test.to_frame().copy()
        df_prevalid['pred'] = Y_pred
                
    #df = pd.Series(df_dict).to_frame().T
    
    if Y_test_surv is not None:
        if evalm == "cox-nloglik":
            Y_pred_binary = np.where(Y_pred < 0.5, 1, 0)
            #hr = get_hr(Y_test_surv, Y_pred_binary)
        #else:   
            #hr = get_hr(Y_test_surv, Y_pred)
            
        #df["HR"] = hr
        
    # return(df, df_prevalid)
    return(df_prevalid)

def unit_train(xgbc, X_train, Y_train, X_test, Y_test, ft_seqs, evalm, HR_calc=False, RFE_criteria="gain", weights=None):
    """Main XGBoost engine that calculates error, vimp, gridsearch, hr
    
    Parameters
    ----------
    xgbc: initialized sklearn.XGBRegressor or sklearn.XGBClassifier 
    X_train: pd.DataFrame
        The input feature dataframe for training
    Y_train: pd.DataFrame
        Label(s) dataframe for training that contains at least one label columns,
        with 2 additonal colunms (survival time and event flag) 
        for calculating hazard ratio
    X_test: pd.DataFrame
        The input feature dataframe for testing
    Y_test: pd.DataFrame
        Label(s) dataframe for testing that contains at least one label columns,
        with 2 optional colunms (survival time and event flag) 
        for calculating hazard ratio 
    ft_seqs: list
        List of numbers for the number of features to used in RFE
    evalm: string
        Evaluation metric used in calculating model score (e.g. auc, accurac, ...)
    HR_calc: bool
        indicate whether to calculate hazard ratio or not
    RFE_criteria: string
        whether to use "gain" or "frequence" for removing features
    weights: np.array
        array of sample weights (default is None)
    
    Returns
    ---------
    df_vimp  : pd.DataFrame
        vimp dataframe (n_features*n_RFE_steps, 
        [feature_name, feature_importance, num_of_feature_used, eval_score])
    df_error : pd.DataFrame
        error metrics (n_RFE,
        [])
    df_grid  : pd.DataFrame
        gridsearch result (n_gridsearch,
        [paramters_used_in_gridsearch ..., f1, size])
    df_prevalid: pd.DataFrame
        contains true label and predicted label (Y_test.shape[1]*n_RFE_steps,
        [true_label, pred, size])
        
    """
    ## initial training
    print(X_train.shape[0], "samples are used.")
    print(X_train.shape[1],"all features are used during initial training.")
    
    ## detect if we need to calculate Hazard Ratio  
    ## disect Y into label and survival time+event_flag
    if HR_calc == "surv":
        Y_train_surv = Y_train[["col_surv", "col_event"]]
        Y_train = Y_train.col_surv_mod
        Y_test_surv = Y_test[["col_surv", "col_event"]]
        Y_test = Y_test.col_surv_mod
        
    if HR_calc == True:
        Y_train_surv = Y_train[["col_surv", "col_event"]]
        Y_train = Y_train.lb
        Y_test_surv = Y_test[["col_surv", "col_event"]]
        Y_test = Y_test.lb
    
    ## check if RFE criteria is properly supplied
    if RFE_criteria not in ["gain", "freq"]:
        print("Your RFE criteria ({}) does not match either 'gain' or 'freq'.".format(RFE_criteria))
        print("using 'gain' as default")
        RFE_criteria = "gain"
    
    ## initialize vimp and prevalidation
    df_vimp = pd.DataFrame()
    df_prevalid = pd.DataFrame()
    df_grid = pd.DataFrame()
    
    features = X_train.columns.values
    ## grid search best hyper-param every 20% reduction in the loop
    for i, k in enumerate(ft_seqs):
        ## Remaining features (sorted)
        top_fts = features[0:k]
        
        if k != len(top_fts):
            raise ValueError("size of selected features does not match RFE_step")
        
        if (i % (round(len(ft_seqs)*0.2)+1) == 0): ## +1 is to avoid error of dividing 0
            print("feature size:", k)
            grid = grid_search(xgbc, X_train.loc[:, top_fts], Y_train, 
                                    X_test.loc[:, top_fts], Y_test, evalm, weights)
            
            ## For admin use only; whether to output just the best grid score or all scores.
            best_grid_score = True  
            ## if True, only output the best grid values. Else, output all search results.
            ## if All search results, the n_estimators value won't be saved unless we also tune this parameter.
            ## GridsearchCV or RandomizedSearchCV only stores the best n_estimators value.
            if best_grid_score:
                best_params = grid.best_params_
                ## extract early_stopping's n_estimators
                best_params["n_estimators"] = grid.best_estimator_.best_iteration
                df_grid_ = pd.DataFrame.from_records([best_params])
                df_grid_["score"] = grid.best_score_
            else:
                ## get all the grid scores (but it does not have n_estimator info for each combination, only the best)
                grid_scores_ = pd.DataFrame(grid.grid_scores_)
                df_grid_ = pd.DataFrame.from_records(grid_scores_.parameters)
                df_grid_["score"] = grid_scores_.mean_validation_score
            
            df_grid_["size"] = k
            df_grid = df_grid.append(df_grid_)
            
        xgbc.set_params(**grid.best_params_)
        xgbc.fit(X_train.loc[:, top_fts], Y_train,
                sample_weight = weights,
                eval_set = [(X_test.loc[:, top_fts], Y_test)], 
                eval_metric = evalm, 
                early_stopping_rounds=100, 
                verbose=False)
        
        if i % 5 == 0:
            evaluation = xgbc.evals_result_['validation_0'][evalm][-1]
            print(evalm+":", evaluation, "at size =", k)
        
        if evalm in ["auc", "merror", "logloss", "mlogloss"]:
            xgbc_cali = CalibratedClassifierCV(xgbc, method='isotonic', cv=5)
            xgbc_cali.fit(X_train.loc[:, top_fts], Y_train, sample_weight = weights)
        else:
            xgbc_cali = None
        
        ## initialize vimp dataframe
        df_vimp_ = pd.DataFrame({
            "feature": top_fts,
            "size": k})
        
        ## frequency and gain
        xgb_obj = xgbc.get_booster()
        vimp_w = xgb_obj.get_score(importance_type="weight")
        vimp_g = xgb_obj.get_score(importance_type="gain")
        
        vimp_w = pd.Series(vimp_w).to_frame(name="weight")
        vimp_w = vimp_w.reset_index() ## move index (feature name) to column
        vimp_w["weight"] = vimp_w["weight"]/vimp_w.weight.sum() ## normalize weight(freq) by dividing the number by 'total number of freqs' ('total number of splits')
        vimp_g = pd.Series(vimp_g).to_frame(name="gain")
        vimp_g = vimp_g.reset_index() ## move index (feature name) to column
        
        df_vimp_wg = vimp_w.merge(vimp_g, on="index", how="outer")
        df_vimp_wg.rename(columns={"index":"feature"}, inplace=True)
        
        ## merge back to df_vimp and get the right column order
        df_vimp_ = df_vimp_.merge(df_vimp_wg, on="feature", how="outer")
        
        ## sort features based on gain
        features = df_vimp_.sort_values('gain', ascending=False).feature.values
        
        ## prevalidation evaluation
        if HR_calc:
            df_prevalid_ = confM(xgbc, X_test.loc[:, top_fts], Y_test, evalm, xgbc_cali, Y_test_surv)
        else:
            df_prevalid_ = confM(xgbc, X_test.loc[:, top_fts], Y_test, evalm, xgbc_cali)
        #df_error["size"] = X_train.shape[1]
        df_prevalid_["size"] = k
        df_prevalid_["n_estimators"] = xgbc.best_iteration
        
        ## append dataframes
        df_vimp = df_vimp.append(df_vimp_)
        df_prevalid = df_prevalid.append(df_prevalid_)
    
    df_vimp = df_vimp[~(df_vimp.weight.isna() | df_vimp.gain.isna())]
            
    return(df_vimp, df_grid, df_prevalid)
