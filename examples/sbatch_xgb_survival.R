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
  project_name = "xgboost_survival_outCV_test",  
  training_data = "data/test_regr_surv_tcga_brca.csv", 
  label_name = NA, 
  sample_ID = "Patient_ID", 
  surv_col = "os_time",
  event_col = "os_event", 
  weight_col = NA, ## colname for the weights
  RFE_step = 10, 
  
  engine = "XGBoost_surv.py", 
  RFE_criteria = 'gain',  ## gain or frequency
  split_CVs = T,
  queue_priority = "short" ## short, medium, long
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
#   mutate(os_time = os_time + 1) ## create a dummy feature
# modified_filename <- "data/df_dummy_modified.csv"

# runSpec$training_data <- modified_filename
# readr::write_csv(df.mod, paste0(runSpec$project_home, "/", runSpec$training_data))

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



