# For R engines: this script can be run in http://rescomp5105.gene.com:8080
# just open this script on RStudio and then click "Source" button 
#

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 2,

  project_home = "~/EVE/tests",
  engine = "lasso.r",
  project_name = "lasso_regression_outCV_test",  
  training_data = "data/test_regr_surv_tcga_brca.csv", 
  label_name = "os_time", #"pam50_RNAseq", 
  sample_ID = "Patient_ID", 
  surv_col = NA, #NA, 
  event_col = NA, #NA, 
  family = "gaussian",
  num_CV = 3, # -1 for LOOCV
  
  ## colname for the weights
  weight_col = NA, 
  
  split_CVs = T,
  queue_priority = "short", ## short, medium, long
  
  ## lass specific parameters
  alpha = 0.5,
  nCv4lambda = 10, #the number of repeated CV to get optimized lambda
  lambdaSum='mean', # after getting nCv4lambda number of lambda, how to summarize them
  lambdaChoice = 'lambda.min' # or lambda.1se, 2 options from cv.glmnet
)
#######################
### Input Data Processing (optional) ###
#######################

## modify data
# In this section, users are supposed to scale numeric features (if necessary); dummify categorical variables ; to provide weight values for each observation.

#Users do not have to remove features with constant values or missing values; they are removed implicitly in the current implementation.

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



