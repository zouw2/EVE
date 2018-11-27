## Instructions ##
## For Python engines: this script must be ran using command line 
# "module load apps/python"
# "Rscript sbatch_xxx.R"
#
# For R engines: this script can be run in http://rescomp5105.gene.com:8080
# just open this script on RStudio and then click "Source" button 
#

################
## User Input ##
################
runSpec <- list(
  seed_base = 1000,
  num_seeds = 2,

  project_home = "~/ml-pipeline/docs/Tutorial",
  engine = "rfeSRCCv3.r", ## ML engine
  #engine = "XGBoost_clasi.py",
  #project_name = "clasi_xgb_outCVs_test4",
  project_name = "survival_rfsrc_outCVs_test4",  
  training_data = "data/TCGA_BCRA_Survival_cleaned_subset.csv", 
  #training_data = "data/TCGA_BCRA_pam50_cleaned_subset.csv",
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
source("~/ml-pipeline/ML_submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



