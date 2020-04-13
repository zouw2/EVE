#module load apps/R
#para = c('1001', "~/2OBD/PDL1mab/go28915_oak/ihc/log/digital_pred_path1/metainfo.Rdata", '0')
#para = c('1002', "~/Biomarker/BDAT/IMpassion130/log/IMpassion130_FMI_PFS_original/metainfo.Rdata", '1')
library(glmnet)
library(survival)
source('~/EVE/eve/models/ml_functions.R')

# if( ! (R.version$major >= 3 && R.version$minor >= 5.1 )) {
#   stop(paste('reporting needs R version 3.5.1'))
# }
# 
# print("EVE commit:")
# git2r::revparse_single(git2r::repository('~/EVE'),"HEAD") # this line works for R 3.4.3/3.5.1 under rescomp3

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

########################
## End of User Inputs ##
########################
if (!file.exists(runSpec$outP)){
  dir.create(runSpec$outP, recursive =T)
} 

if(runSpec$family == "cox"){
  runSpec$label_name <- runSpec$surv_col
}

outF <- paste0(runSpec$outP, 
               "/glmnet_rdata_", runSpec$label_name, 
               "_seed", specLocal$seed, 
               "_cv", specLocal$cv_id_curr, ".rdata")
#outF <- paste(runSpec$outP, '/glmnet_', paste(unlist(specLocal), collapse="_"),'.rdata', sep='')

if( file.exists(outF) && T) {
  print(paste(outF,'already exits. quitting...'))
  q()
}  

#########################
## read & process data ##
#########################
if(grepl('csv$', tolower(runSpec$training_data))) df <- read.csv(paste(runSpec$project_home, runSpec$training_data, sep='/'))

if(grepl('rds$', tolower(runSpec$training_data))) df <- as.data.frame(readRDS(paste(runSpec$project_home, runSpec$training_data, sep='/')))

if(is.na(runSpec$sample_ID)|is.null(runSpec$sample_ID)|(runSpec$sample_ID=="")){
  print("Use index as sample ID")
  runSpec$sample_ID <- "RowIndex"
} else {
  rownames(df) <- as.character(df[[runSpec$sample_ID]])
  stopifnot(runSpec$sample_ID %in% colnames(df))
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
if (runSpec$family %in% c("multinomial", "binomial")) {
  stopifnot(runSpec$label_name %in% colnames(df))
  
  df[[runSpec$label_name]] <- factor(df[[runSpec$label_name]])
  y <- df[, runSpec$label_name, drop=F] # y will be handled similarly for cox and categorical outcome
} 
if (runSpec$family %in% c("gaussian")) {
  stopifnot(runSpec$label_name %in% colnames(df))
  y <- df[, runSpec$label_name, drop=F] # y will be handled similarly for cox and categorical outcome
}  

if(any(is.na(as.vector(y)))){
  print(paste('there are', sum(is.na(as.vector(y))),'missing values in the outcome data file'))
}

col2drop <- c(runSpec$label_name, runSpec$sample_ID,
              runSpec$surv_col,   runSpec$event_col, runSpec$weight_col, runSpec$variable2ignore)

x <- data.matrix( df[,!(colnames(df) %in% col2drop)] ) ## why data.matrix not data.frame

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

# not applicable for glmnet

###############
## handle CV ##
###############
RNGversion('3.5.1') # backward compatible with earlier eve runs in rescomp3
set.seed(specLocal$seed)

if (runSpec$family %in% 'cox') {
  cvList <- caret::createFolds(y[, "col_surv"] * (as.integer(y[, "col_event"]) - 0.5), 
                               k = runSpec$num_CV, 
                               list = TRUE, returnTrain = FALSE)
}

if (runSpec$family %in% c("multinomial", "binomial")) {
  stopifnot(is.factor(y[, runSpec$label_name]) )
  print(table(y[, runSpec$label_name]))
  cvList <- caret::createFolds(y[, runSpec$label_name], 
                               k = runSpec$num_CV, 
                               list = TRUE, returnTrain = FALSE)
}

if (runSpec$family %in% c("gaussian")) {
  stopifnot(is.numeric(y[, runSpec$label_name]))
  cvList <- caret::createFolds(y[, runSpec$label_name], 
                               k = runSpec$num_CV, 
                               list = TRUE, returnTrain = FALSE)
}

## Note: using for loop is because we can keep track of cv_id in the output

df_vimp <- data.frame()
df_pred <- data.frame()

cv_id <- 1

start.time <- Sys.time()

for (cv.idx in cvList){
  per_cv_seed <- (specLocal$seed - runSpec$seed_base) * (nrow(x) + 5) + cv_id #this seed should be different per cv per split/seed
  
  if (specLocal$cv_id_curr == 0){
    print("CVs within job")
    print(paste("Fold number:", cv_id))
    X_train = x[-cv.idx, featureList, drop=F]
    Y_train = y[-cv.idx, , drop=F]
    X_test  = x[ cv.idx, featureList, drop=F]
    Y_test  = y[ cv.idx, , drop=F]
    #Y_test  = cbind(Y_test, rownames(Y_test))
    #colnames(Y_test) <- c(colnames(y), runSpec$sample_ID)
    #Y_test[runSpec$sample_ID] <- rownames(Y_test)
    print(paste('using', paste(head(cv.idx, 10), collapse=','),',etc, as validation' ))
    
    df.out <- glmnetCVwrapper2(X_train, Y_train, X_test, Y_test, 
                               seed=per_cv_seed, 
                               glmnetFam = runSpec$family, 
                               a1 = runSpec$alpha, 
                               nCv4lambda = runSpec$nCv4lambda, 
                               lambdaSum = match.fun(runSpec$lambdaSum), runPairs=runSpec$runPairs,
                               lambdaChoice = runSpec$lambdaChoice, 
                               w = runSpec$weight.value[-cv.idx])
    if(nrow(df.out$features) > 0){
      df_vimp_tmp <- data.frame( df.out$features,
                                "lambda" = unname(df.out$lambda),'alpha' =  df.out$alpha,
                                "cv" = cv_id, stringsAsFactors = F)
      colnames(df_vimp_tmp) <- c("feature", 'coef', ifelse(ncol(df.out$features) == 2, NULL, 'vimp'), "lambda", 'alpha','cv')
      
      df_vimp <- rbind(df_vimp, df_vimp_tmp)
    }
    df_pred_tmp <- df.out$pred
    df_pred_tmp["cv"] <- cv_id
    df_pred_tmp[[runSpec$sample_ID]] <- rownames(X_test) ## add sample ID here
    
    
    df_pred <- rbind(df_pred, df_pred_tmp)
    
  } else {
    if (specLocal$cv_id_curr == cv_id){
      print("1 CV per job")
      print(paste("Fold number:", cv_id))
      X_train = x[-cv.idx, featureList, drop=F]
      Y_train = y[-cv.idx, , drop=F]
      X_test  = x[ cv.idx, featureList, drop=F]
      Y_test  = y[ cv.idx, , drop=F]
      #Y_test  = cbind(Y_test, rownames(Y_test))
      #colnames(Y_test) <- c(colnames(y), runSpec$sample_ID)
      #Y_test[runSpec$sample_ID] <- rownames(Y_test)
      print(paste('using', paste(head(cv.idx, 10), collapse=','),',etc, as validation' ))
      df.out <- glmnetCVwrapper2(X_train, Y_train, X_test, Y_test, 
                                 seed=per_cv_seed, 
                                 glmnetFam = runSpec$family, 
                                 a1 = runSpec$alpha, 
                                 nCv4lambda = runSpec$nCv4lambda, 
                                 lambdaSum = match.fun(runSpec$lambdaSum), runPairs=runSpec$runPairs,
                                 lambdaChoice = runSpec$lambdaChoice, 
                                 w = runSpec$weight.value[-cv.idx])
      
      if(nrow(df.out$features) > 0){
        df_vimp_tmp <- data.frame( df.out$features,
                                "lambda" = df.out$lambda, 'alpha' =  df.out$alpha,
                                "cv" = cv_id, stringsAsFactors = F)
        df_vimp <- rbind(df_vimp, df_vimp_tmp)
      }else{
        cat('no feature retained for split', cv_id,'under local seed', specLocal$seed)
      }
      
      df_pred_tmp <- df.out$pred
      df_pred_tmp["cv"] <- cv_id
  #    df_pred_tmp[[runSpec$sample_ID]] <- rownames(X_test) ## add sample ID here
      df_pred_tmp[[runSpec$sample_ID]] <- rownames(df_pred_tmp)
      rownames(df_pred_tmp) <- NULL
      
      df_pred <- rbind(df_pred, df_pred_tmp)
    }
  }
  cv_id <- cv_id + 1
}

## added 'size' column to work with reporting functions
## note: plotting performance against # of feature retained is not meaningful for lasso.
## note: size for lasso is just the total # of features in the input matrix

if(nrow(df_vimp) > 0) df_vimp$size <- dim(x)[2]
df_pred$size <- dim(x)[2]

stopifnot(!any(duplicated(df_pred$sample_ID)))
save(df_pred, df_vimp, file=outF)

end.time <- Sys.time()
time.taken <- end.time - start.time

print(paste('time taken: ', time.taken)) # give the running length helps user to pick a queue type