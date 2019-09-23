## To submit the job, simple run this script, 
## or click "Source" button if you are using RStudio

#modified from sbatch_lasso_binary.R

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 4,
  num_CV = 5, # -1 for LOOCV

  project_home = "~/EVE/examples",
  project_name = "lasso_binary2",  
  training_data = "data/test_binaryclass_tcga_brca.csv", # or use saveRDS to save rds file
  label_name = "pam50_RNAseq",
  sample_ID = "Patient_ID", 
  surv_col = NA, 
  event_col = NA, 
  weight_col = NA, # colname for the weights 
  
  engine = "lasso.r",
  family = "binomial",
  split_CVs = T,
  queue_priority = "short", ## short, medium, long
  
  ## lasso specific parameters
  alpha = c(0, 0.5, 1),                 # use 0.5 for elastic net, or 0 for ridge
  nCv4lambda = 10,           # the number of repeated CV to get optimized lambda
  lambdaSum='mean', # after getting nCv4lambda number of lambda, how to summarize them
  lambdaChoice = 'lambda.min' # or lambda.1se, 2 options from cv.glmnet
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
#   mutate(pam50_RNAseq2 = pam50_RNAseq + 1) ## create a dummy feature
# modified_filename <- "data/df_dummy_modified.csv"

# runSpec$training_data <- modified_filename
# readr::write_csv(df.mod, paste0(runSpec$project_home, "/", runSpec$training_data))

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



