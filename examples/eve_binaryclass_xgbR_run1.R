## To submit the job, simple run this script, 
## or click "Source" button if you are using RStudio

################
## User Input ##
################
runSpec <- list( # to be moved to  sbatch later
  seed_base = 1000,
  num_seeds = 4,
  num_CV = 3, # -1 for LOOCV
  
  project_home = "~/EVE/examples",
  project_name = "xgboostR_binary_1",  
  training_data = "data/test_binaryclass_tcga_brca.csv", 
  label_name = "pam50_RNAseq", 
  sample_ID = "Patient_ID", 
  surv_col = NA, 
  event_col = NA, 
  weight_col = NA, ## colname for the weights
  nCv4lambda = 3,  # number of CV used for xgb.cv()
  RFE_step = 20.2 , 
  
  engine = "XGBoost.r", 
  family = "binomial",
  RFE_criteria = 'gain',
  split_CVs = T,
  queue_priority = "short", ## short, medium, long
#  server_management =  '-n 12 -R "span[hosts=1] rusage[mem=12]"', # this asks for a cluster with 12 o more CPU, each with 12G memory or more
nthread=16, # this will feed into R implementation xgb.cv and xgb.train. a large number may increase the speed of cv, but make you wait longer to get access to the requested node. the python implementation of xgboost seems to able to grab what ever nodes avaiable; for the R version xgboost availabe in R 3.5.1, it seems we have to explicitly request it.
  "n_estimators" = 3000, # the largest number of round
#  'early_stopping_rounds' = 200, # how often to check early stopping, if not specified, use 1/eta
  evalm = NULL, # this info will determine the eval_metric element in xgb_params. if not specified, the default will be ('binomial' = "logloss", "multinomial"= "mlogloss", "gaussian" = 'rmse'). customer evaluation functions currently implement: f1_harmonic2,f1_meanlog2 (similar to geometric mean),f1_arithmetic2
  tune_params = list(eta = c(1e-3, 1e-2, 0.1), gamma =c(0, 10, 20), max_depth = c(1, 6, 10), max_delta_step=c(0,3,5),colsample_bytree =c(1, 0.7, 0.5) ,lambda =c(0, 1, 10), alpha=c(0, 1,10) ), # if a vector of length 1 is specified, the parameter will be assigned to xgb (and potential overwrite default parameter if conflict)
  pct_feature_for_grid_search = seq(1, 0.2, -0.3), # this is c(1, 0.7, 0.4), parameter tuning with xgb.cv() will be repeated when approximately 100%, 70% and 40% of total features are kept 
  max_num_grid = 30 # if there are more than this number of combinations, pick a random portion of it. too many grid search will make code very slow
  
)


########################################
### Input Data Processing (optional) ###
########################################

## modify data
# XGBoost engine runs in python which does not understand characters,
# therefore, user has to convert all categorical/characters into numeric. 
# Possibly dummify categorical/cardinal variables
# (sometimes encode each cardinal values into numerical number works just fine, but no guarantee).
# Provide weight values for each observation in a column if necessary.

# Users do not have to remove features with constant values or missing values; 
# XGboost has its own method of encoding missing values and remove low variance feature.

## Example of how to modify input data
# library(dplyr)
# df <- readr::read_csv(paste0(runSpec$project_home, "/", runSpec$training_data))
# df.mod <- df %>% 
#   mutate(pam50_RNAseq = pam50_RNAseq + 1) ## create a dummy feature
# modified_filename <- "data/df_dummy_modified.csv"

# runSpec$training_data <- modified_filename
# readr::write_csv(df.mod, paste0(runSpec$project_home, "/", runSpec$training_data))

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



