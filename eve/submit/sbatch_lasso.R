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
  engine = "lasso1.r", ## ML engine
  project_name = "survival_lasso_inCVs_test2",  
  training_data = "data/TCGA_BCRA_Survival_cleaned_subset.csv", 
  #training_data = "data/TCGA_BCRA_pam50_cleaned_subset.csv",
  label_name = NA, #"pam50_RNAseq", 
  sample_ID = "Patient_ID", 
  surv_col = "os_time", #NA, 
  event_col = "os_event", #NA, 
  family  = 'cox', # this string will be fed to glmnet(family = ). "cox", "multinomial", "binomial"
  #referenece_level  = 'DESERT',
  num_CV = 3, # -1 for LOOCV
  RFE_step = 10, 
  
  ## change the weight of positive label
  weight_col = NA, 
  
  split_CVs = F,
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

#######################
## End of User Input ##
#######################
source("~/EVE/eve/submit/utils/sbatch_submit.R")
sbatch_submit(runSpec)



