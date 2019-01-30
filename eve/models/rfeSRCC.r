#module load apps/R
#para = c('1001', "~/ml-pipeline/docs/Tutorial/log/survival_rfsrc_outCVs_test4/metainfo.Rdata", '0')
source('~/EVE/eve/models/ml_functions.R')

#################
## User Inputs ##
################# 
para <- commandArgs(trailingOnly = TRUE)
print(para)

load(para[2]) ## get runSpec
seed <- as.integer(para[1])
cv_id_curr <- as.integer(para[3])

########################
## End of User Inputs ##
########################
if (!file.exists(runSpec$outP)){
  dir.create(runSpec$outP, recursive =T)
} 

## check whether results files already exist,
## if so, don't run 
## this may cause problem if we want to 'overwrite' the output
# specLocal <- list(
#   seed = seed,
#   cv_id_curr = cv_id_curr
# )
# 
# outF <- paste0(runSpec$outP, '/rfsrc_', paste(unlist(specLocal), collapse="_"), '.rdata')
# if( file.exists(outF) && T) {
#   print(paste(outF,'already exits. quitting...'))
#   q()
# }

#########################
## read & process data ##
#########################
df <- read.csv(paste(runSpec$project_home, runSpec$training_data, sep='/') )

if(is.na(runSpec$sample_ID)|is.null(runSpec$sample_ID)|(runSpec$sample_ID=="")){
  print("Use index as sample ID")
  runSpec$sample_ID <- "RowIndex"
} else {
  rownames(df) <- as.character(df[[runSpec$sample_ID]])
  stopifnot(runSpec$sample_ID %in% colnames(df))
  stopifnot(!any(duplicated(as.character(df[[runSpec$sample_ID]]))))
}

col2drop <- c(runSpec$label_name, runSpec$sample_ID,
              runSpec$surv_col,   runSpec$event_col, runSpec$weight_col)

x <- df[,!(colnames(df) %in% col2drop)]
y <- df[,c(runSpec$surv_col, runSpec$event_col)]
colnames(y) <- c("col_surv", "col_event")

###########################
## handle input features ##
###########################

## ToDo: allow user to input feature of interest
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
} else if(length(RFE_step) == 1){ 
  stopifnot(runSpec$RFE_step <= nFeature)
  
  ## use all features
  if(RFE_step == 0){ 
    sizes <- nFeature 
  ## use fraction of features
  } else if (0.0 < RFE_step & RFE_step < 1){
    
    ft <- nFeature
    sizes <- c()
    step_size <- 1
    
    while (step_size >= 1){
      step_size <- round(RFE_step*ft)
      sizes <- c(sizes, ft) 
      ft <- ft - step_size
    }
  }
  ## use fixed step_size
  else {
    sizes <- seq(nFeature, 1,  runSpec$RFE_step*-1)
  }
}
print("RFE steps:")
print(sizes)

###############
## handel CV ##
###############
set.seed(seed)
cvList <- caret::createFolds(y$col_surv * (as.integer(y$col_event) - 0.5), 
                             k = runSpec$num_CV, 
                             list = TRUE, returnTrain = FALSE)

df_vimp <- data.frame()
df_pred <- data.frame()
cv_id <- 1

start.time <- Sys.time()

for (cv.idx in cvList){
  if (cv_id_curr == 0){
    print("CVs within job")
    print(paste("Fold number:", cv_id))
    X_train = x[-cv.idx, featureList, drop=F]
    Y_train = y[-cv.idx, ]
    X_test  = x[ cv.idx, featureList, drop=F]
    Y_test  = y[ cv.idx, ]
    print(paste('using', paste(head(cv.idx, 10), collapse=','),',etc, as validation' ))
    df.out <- rfeSRCCv3(X_train, Y_train, X_test, Y_test, sizes, seed)
    
    df_vimp_tmp <- df.out$df_vimp
    df_pred_tmp <- df.out$df_pred
    
    df_vimp_tmp["cv"] <- cv_id
    df_pred_tmp["cv"] <- cv_id
    df_pred_tmp[[runSpec$sample_ID]] <- rownames(X_test)
    print(paste0("dim df_vimp: ", dim(df_vimp_tmp)))
    print(paste0("dim df_pred: ", dim(df_pred_tmp)))
    df_vimp <- rbind(df_vimp, df_vimp_tmp)
    df_pred <- rbind(df_pred, df_pred_tmp)
    
  } else if (cv_id_curr == cv_id){
    print("1 CV per job")
    print(paste("Fold number:", cv_id))
    X_train = x[-cv.idx, featureList, drop=F]
    Y_train = y[-cv.idx, ]
    X_test  = x[ cv.idx, featureList, drop=F]
    Y_test  = y[ cv.idx, ]
    print(paste('using', paste(head(cv.idx, 10), collapse=','),',etc, as validation' ))
    df.out <- rfeSRCCv3(X_train, Y_train, X_test, Y_test, sizes, seed)
    
    df_vimp_tmp <- df.out$df_vimp
    df_pred_tmp <- df.out$df_pred
    
    df_vimp_tmp["cv"] <- cv_id
    df_pred_tmp["cv"] <- cv_id
    df_pred_tmp[[runSpec$sample_ID]] <- rownames(X_test)
    
    df_vimp <- rbind(df_vimp, df_vimp_tmp)
    df_pred <- rbind(df_pred, df_pred_tmp)
  }
  cv_id <- cv_id + 1
}

################
## save files ##
################
## if results directory not exist, create one
if(!dir.exists(runSpec$outP)){
  dir.create(runSpec$outP, recursive=TRUE)
}

#out_vimp   <- paste0(runSpec$outP, "/rfsrc_vimp_", runSpec$surv_col, "_seed", seed, "_cv", cv_id_curr, ".csv")
#out_preval <- paste0(runSpec$outP, "/rfsrc_prevalidation_", runSpec$surv_col, "_seed", seed, "_cv", cv_id_curr, ".csv")
#write.csv(df_vimp, out_vimp, row.names = F)
#write.csv(df_pred, out_preval, row.names = F)

outF <- paste0(runSpec$outP, "/rfsrc_rdata_", runSpec$surv_col, "_seed", seed, "_cv", cv_id_curr, ".rdata")
save(df_pred, df_vimp, file=outF)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste('time taken:', time.taken)) # give the running length helps user to pick a queue type