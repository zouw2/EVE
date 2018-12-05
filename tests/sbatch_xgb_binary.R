## Instructions ##
## For Python engines: this script must be ran using command line 
# "module load apps/python"
# "Rscript ~/EVE/tests/sbatch_xgb_binary.R"

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 2,

  project_home = "~/EVE/tests",
  engine = "XGBoost_clasi.py", 
  project_name = "xgboost_binary_outCV_test",  
  training_data = "data/test_binaryclass_tcga_brca.csv", 
  label_name = "pam50_RNAseq", 
  sample_ID = "Patient_ID", 
  surv_col = NA, 
  event_col = NA, 
  num_CV = 3, # -1 for LOOCV
  RFE_step = 10, 
  
  ## change the weight of positive label
  scale_pos_weight = 0.3, 
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



