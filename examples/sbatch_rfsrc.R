## To submit the job, simple run this script, 
## or click "Source" button if you are using RStudio

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 2,
  num_CV = 3, # -1 for LOOCV

  project_home = "~/EVE/tests",
  project_name = "rfsrc_outCV_test",  
  training_data = "data/test_regr_surv_tcga_brca.csv", 
  label_name = NA, 
  sample_ID = "Patient_ID", 
  surv_col = "os_time",  
  event_col = "os_event", 
  weight_col = NA, ## not yet implemented for rfsrc engine
  RFE_step = 10, 
  
  engine = "rfeSRCC.r", ## ML engine
  RFE_criteria = 'permute',  
  split_CVs = T,
  queue_priority = "short" ## short, medium, long
)
########################################
### Input Data Processing (optional) ###
########################################

## modify data
# In this section, users are supposed to dummify categorical variables ; 
# Remaining character variables will be converted to factors and then integers, which is usually not what we want.
# User also has to ensure there is no NA in the data.

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



