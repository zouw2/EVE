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
  project_name = "xgboost_binary_outCV_test_v2",  
  training_data = "data/test_binaryclass_tcga_brca_weightcol.csv", 
  label_name = "pam50_RNAseq", 
  sample_ID = "Patient_ID", 
  surv_col = NA, 
  event_col = NA, 
  num_CV = 3, # -1 for LOOCV
  RFE_step = 10, 
  
  ## change the weight of positive label
  #scale_pos_weight = 0.3, 
  weight_col = "weight_col",
  RFE_criteria = 'gain',  
  
  split_CVs = T,
  queue_priority = "short" ## short, medium, long
  
)
#######################
### Input Data Processing (optional) ###
#######################

#df <- read.csv(paste0(runSpec$project_home, "/", "data/test_binaryclass_tcga_brca.csv"))
#df$weight_col <- ifelse(df$pam50_RNAseq==1, 0.3, 1)
#write.csv(df, paste0(runSpec$project_home, "/", runSpec$training_data), 
#          row.names = F)

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



