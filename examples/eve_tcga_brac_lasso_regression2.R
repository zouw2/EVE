## To submit the job, simple run this script, 
## or click "Source" button if you are using RStudio

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 4,
  num_CV = 3, # -1 for LOOCV

  project_home = "~/EVE/examples",
  project_name = "lasso_regression_noSplitCV",  
  training_data = "data/test_regr_surv_tcga_brca.csv", # or use saveRDS to save rds file
  label_name = "os_time", 
  sample_ID = "Patient_ID", 
  surv_col = NA, 
  event_col = NA,  
  weight_col = NA, # colname for the weights
  
  engine = "lasso.r",
  family = "gaussian",
  split_CVs = F,
  queue_priority = "short", ## short, medium, long
  
  ## lasso specific parameters
  alpha = 1,                 # use 0.5 for elastic net, or 0 for ridge
  nCv4lambda = 10,           # the number of repeated CV to get optimized lambda
  lambdaSum='mean', # after getting nCv4lambda number of lambda, how to summarize them
  lambdaChoice = 'lambda.min', # or lambda.1se, 2 options from cv.glmnet
  variable2ignore =c()
)
########################################
### Input Data Processing (optional) ###
########################################

## modify data
# In this section, users are supposed to scale numeric features (if necessary); 
# dummify categorical variables ; to provide weight values for each observation. 
# Remaining character variables will be converted to factors and then integers, which is usually not what we want.

# Users do not have to remove features with constant values or missing values; 
# they are removed implicitly in the current implementation.

## Example of how to modify input data
# library(dplyr)
# df <- readr::read_csv(paste0(runSpec$project_home, "/", runSpec$training_data))
# df.mod <- df %>% 
#   mutate(os_time = os_time + 1) ## create a dummy feature
# modified_filename <- "data/df_dummy_modified.csv"

# runSpec$training_data <- modified_filename
# readr::write_csv(df.mod, paste0(runSpec$project_home, "/", runSpec$training_data))

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



