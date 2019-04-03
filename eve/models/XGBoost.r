#module load apps/R # xgboost in R3.4.3 is fairly old version: "0.6.4.1"
#module load apps/R/3.5.1-Bioc-3.8-test/test

#para = c('1002', "~/2OBD/PDL1mab/go28915_oak/clinical_unofficial/log//metainfo.Rdata", '1')



source('~/EVE/eve/models/ml_functions.R')


#################
## User Inputs ##
################# 
para <- commandArgs(trailingOnly = TRUE)
print(para)

load(para[2]) ## get runSpec

specLocal <- list(
  seed = as.integer(para[1]), 
  cv_id_curr = as.integer(para[3])
)

runSpec[['curr_cv']] <- specLocal$cv_id_curr
runSpec[['localSeed']] <- specLocal$seed

# outP is added to runSpec in sbatch_submit() 
#runSpec$outP <- paste0(runSpec$project_home, "/results/", runSpec$project_name)

########################
## End of User Inputs ##
########################
if (!file.exists(runSpec$outP)){ 
  dir.create(runSpec$outP, recursive =T)
} 

# if(runSpec$family == "cox"){
#   runSpec$label_name <- runSpec$surv_col
# }

outF <- paste0(runSpec$outP, 
               "/xgbR_", runSpec$label_name, 
               "_seed", specLocal$seed, 
               "_cv", specLocal$cv_id_curr, ".rdata")

if( file.exists(outF) && T) {
  print(paste(outF,'already exits. quitting...'))
  q()
}  

#########################
## read & process data ##
#########################
if(grepl('csv$', tolower(runSpec$training_data))) df <- read.csv(paste(runSpec$project_home, runSpec$training_data, sep='/'))

if(grepl('rds$', tolower(runSpec$training_data))) df <- readRDS(paste(runSpec$project_home, runSpec$training_data, sep='/'))

if(is.na(runSpec$sample_ID)|is.null(runSpec$sample_ID)|(runSpec$sample_ID=="")){
  print("Use index as sample ID")
  runSpec$sample_ID <- "RowIndex"
} else {
  stopifnot(runSpec$sample_ID %in% colnames(df))
  rownames(df) <- as.character(df[[runSpec$sample_ID]])
  stopifnot(!any(duplicated(as.character(df[[runSpec$sample_ID]]))))
}

if( is.null(runSpec$weight_col) || is.na(runSpec$weight_col) || runSpec$weight_col==""){
  runSpec$weight.value <- rep(1, nrow(df))
}else{
  stopifnot(runSpec$weight_col %in% colnames(df))
  runSpec$weight.value <- df[, runSpec$weight_col]
}

if (runSpec$family %in%  'cox') {
  stopifnot(all(c(  runSpec$surv_col, runSpec$event_col) %in% colnames(df)))
  y <- data.matrix( df[,c(runSpec$surv_col, runSpec$event_col)] )
  colnames(y) <- c("col_surv", "col_event") 
  
}  

#The XGBoost algorithm requires that the class labels (Site names) start at 0 and increase sequentially
#https://rpubs.com/mharris/multiclass_xgboost

if (runSpec$family %in% c("multinomial", "binomial")) {
  stopifnot(runSpec$label_name %in% colnames(df))
  
  df[[runSpec$label_name]] <- as.integer(factor(df[[runSpec$label_name]])) - 1
  y <- df[, runSpec$label_name, drop=F] # y will be handled similarly for cox and categorical outcome
} 

if (runSpec$family %in% c("gaussian")) {
  stopifnot(runSpec$label_name %in% colnames(df))
  y <- df[, runSpec$label_name, drop=F] # y will be handled similarly for cox and categorical outcome
}  

# the reason to keep y is to handle survival outcome (which has 2 columns). So far xgboost in R does not handle survival data yet; and its documuation wants label as vector. 
 

if(any(is.na(as.vector(y)))){
  print(paste('there are', sum(is.na(as.vector(y))),'missing values in the outcome data file'))
}


col2drop <- c(runSpec$label_name, runSpec$sample_ID,
              runSpec$surv_col,   runSpec$event_col, runSpec$weight_col)

x <- data.matrix( df[,!(colnames(df) %in% col2drop)] ) 

###########################
## handle input features ##
###########################


featureList <- NULL
if(is.null(featureList)){
  featureList <- colnames(x)
}else{
  stopifnot(all(featureList %in% colnames(x)))
}

nFeature <- length(featureList)
######################
## handle RFE steps ##
######################

RFE_step <- runSpec$RFE_step

## if user provides a vector of steps
if(length(RFE_step) > 1){ 
  sizes <- RFE_step
  
  ## if user provides a single value
} else {
  sizes <- decideRFEseq(RFE_step = runSpec$RFE_step, ft =nFeature)
} 


stopifnot(max(sizes) <= nFeature)

print("RFE steps:")
print(sizes)

###############
## parameters ##
###############

# default parameters
# a good source of xgboost documentation: https://sites.google.com/view/lauraepp/parameters

xgb_params <- list(
  
#  "n_estimators" = 3000, # this is not a parameter for R
  "eta" = 0.01, # 0.3, [0,1], usually 0.01 - 0.2
  #  "gamma" = 0, # 0, [0,∞], minimum loss reduction required to make a further partition on a leaf node of the tree. The larger, the more conservative the algorithm will be.
  #  "max_depth" = 4, # 6, 0 means no limit. [0,∞]
  #  "min_child_weight" = 1, # 1, [0,∞], If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. In linear regression mode, this simply corresponds to minimum number of instances needed to be in each node. The larger, the more conservative the algorithm will be.
  #  "max_delta_step"= 0, # 0, [0,∞], Usually this parameter is not needed, but it might help in logistic regression when class is extremely imbalanced. Set it to value of 1-10 might help control the update
  #  "subsample"= 1, # 1, [0,1], subsample ratio of the training instance. Setting it to 0.5 means that XGBoost randomly collected half of the data instances to grow trees and this will prevent overfitting.
  "colsample_bytree"= 0.8, # 1, [0,1], subsample ratio of columns when constructing each tree
  
  # although the r document says lambda and alpha are parameters only for linear booster in r; the general xgboost document says otherwise. From this discussion (https://github.com/dmlc/xgboost/issues/2656), both alpha and lambda are usable for gbtree. 
  #  "lambda"= 1, # 1 (python) or 0 (r), L2 regularization term on weights, increase this value will make model more conservative
  #  "alpha"= 0, # 0, L1 regularization term on weights, increase this value will make model more conservative
  eval_metric = ifelse ( !(is.null( runSpec$evalm ) || is.na(runSpec$evalm)), runSpec$evalm, unname(c('binomial' = "logloss", "multinomial"= "mlogloss", "gaussian" = 'rmse')[runSpec$family])) ,
  #  "silent": 1,
  objective = unname( c('binomial' = "binary:logistic", "multinomial"= "multi:softprob", "gaussian" = 'reg:linear')[runSpec$family] )  
 
  #  "njobs": -1,
  #  "random_state": seed
)


if(length(unique(y[, 1])) > 2) {
  xgb_params[["num_class"]] = length(unique(y[, 1]))
  stopifnot(runSpec$family == "multinomial")
}  

###############
## handle CV ##
###############
set.seed(specLocal$seed)

# the following line will be a problem for survival outcome
cvList <- stratFold(y[, 1], runSpec$family,  runSpec$num_CV)


stopifnot( specLocal$cv_id_curr %in% 1:length(cvList))
cv.idx <- cvList[[specLocal$cv_id_curr]]
 
    print("1 CV per job")

    X_train = x[-cv.idx, featureList, drop=F]
    
# the reason to keep y is to handle survival outcome (which has 2 columns). So far xgboost in R does not handle survival data yet; and its documuation wants label as vector. 
    Y_train = y[-cv.idx, 1]
    X_test  = x[ cv.idx, featureList, drop=F]
    Y_test  = y[ cv.idx, 1]
 
    print(paste('using', paste(head(cv.idx, 10), collapse=','),',etc, as validation' ))
    df.out <- xgbCVwrapper(X_train, Y_train, X_test, Y_test, ft_seqs = sizes, 
                          # glmnetFam = runSpec$family, 
                          # nCv4lambda = runSpec$nCv4lambda, 
                          weight_train = runSpec$weight.value[-cv.idx], params = xgb_params, spec = runSpec)
    
   df.out <- lapply(df.out, function(x) {
     row.names(x) <- NULL
     x})
   df_pred =df.out[[3]]; df_vimp = df.out[[1]]; df_grid= df.out[[2]]
   
   
  save(df_pred, df_vimp, df_grid, file=outF)
  
  sessionInfo()