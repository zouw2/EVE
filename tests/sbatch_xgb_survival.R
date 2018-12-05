## Instructions ##
## For Python engines: this script must be ran using command line 
# "module load apps/python"
# "Rscript ~/EVE/tests/sbatch_xgb_survival.R"

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 2,

  project_home = "~/EVE/tests",
  engine = "XGBoost_surv.py", ## ML engine
  project_name = "xgboost_survival_outCV_test",  
  training_data = "data/test_regr_surv_tcga_brca.csv", 
  label_name = NA, #"pam50_RNAseq", 
  sample_ID = "Patient_ID", 
  surv_col = "os_time", #NA, 
  event_col = "os_event", #NA, 
  num_CV = 3, # -1 for LOOCV
  RFE_step = 10, 
  
  ## change the weight of positive label
  scale_pos_weight = 1, 
  RFE_criteria = 'gain',  
  
  split_CVs = T,
  queue_priority = "short" ## short, medium, long
  
)
#######################
### Input Data Processing (optional) ###
#######################

## modify data

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



