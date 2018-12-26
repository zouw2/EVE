findConstantVar <- function(ds){
  stopifnot(is.matrix(ds) || is.data.frame(ds))
  
  monov <- sapply(colnames(ds), function(v) {
    x <- ds[, v]
    if(length(unique(x)) == 1){
      print(paste('dropping',v,'as it has only 1 possible value'))
      return(v)
    }
    
    if(any(is.na(x))){
      print(paste('dropping', v, 'as it has',sum(is.na(x)),'missing values'))
      return(v)
    }
    
    NA
  })
  monov[!is.na(monov)]
}  

glmnetCVwrapper2 <- function(X_train , Y_train, X_test, Y_test,
                             seed = 27519, 
                             glmnetFam="binomial", a1=1, 
                             nCv4lambda=10, lambdaSum=mean, 
                             lambdaChoice = 'lambda.min',  w=1, ... ){
  
  set.seed(seed)
  stopifnot(is.matrix(X_train)); stopifnot(is.matrix(X_test))
  stopifnot(rownames(X_test) == rownames(Y_test))
  
  if (glmnetFam %in% c("multinomial", "binomial") ) {
    if (class(Y_train) %in% c('matrix','data.frame')) Y_train <- Y_train[, 1]
    #   if (class(Y_test) %in% c('matrix','data.frame'))  Y_test <-  Y_test[, 1]
    stopifnot( is.factor(Y_train) );
    #    stopifnot( is.factor(Y_test) );
    print('training data label counts')
    print(table(Y_train))
  }
  
  if (glmnetFam %in% c("gaussian") ) {
    if (class(Y_train) %in% c('matrix','data.frame')) Y_train <- Y_train[, 1]
    stopifnot( is.numeric(Y_train) );
    print('training data outcome')
    print(summary(Y_train))
  }
  
  if (glmnetFam %in% "cox") {
    
    stopifnot(is.matrix(Y_train) && ncol(Y_train) == 2 && all ( Y_train[, 2] %in% c(T,F)))
    stopifnot(is.matrix(Y_test) && ncol(Y_test) == 2 && all ( Y_test[, 2] %in% c(T,F)))
    colnames(Y_train) <- c('time','status')
    colnames(Y_test) <- c('time','status')
    
    ## drop samples if their time is 0
    idx2drop <- Y_train[, "time"] == 0
    if(any(idx2drop)){
      print(paste("Drop", sum(idx2drop), "samples with survival time = 0."))
      Y_train <- Y_train[!idx2drop, ]
      X_train <- X_train[!idx2drop, ]
    }
  }
  
  if (is.null(w)){ ## if no weight vector is given
    w = rep(1, dim(X_train)[1])
  }
  print(paste('alpha=', a1))
  
  print('dim of training x')
  print(dim(X_train))
  
  monov <- findConstantVar(X_train)
  
  x1 <- X_train[, setdiff(colnames(X_train), monov)]

  maxL <- lambdaSum( sapply(1:nCv4lambda, function(i) {
    set.seed(seed +100 + i)
    
    r1 <- cv.glmnet(x=x1, y=(function(yobj){
      if(glmnetFam =='cox') return(Surv(yobj)); 
      return(yobj) })(Y_train), family=glmnetFam, alpha = a1 , standardize =T, weights=w, ...)
    
    r1[[lambdaChoice]]
    
  }) , na.rm = T)

  g1 <- glmnet(x=x1, y = (function(yobj){
    if(glmnetFam =='cox') return(Surv(yobj)); 
    return(yobj) })(Y_train),  family = glmnetFam, alpha = a1, lambda=maxL, standardize =T, weights=w, ...)
  
  pred.response <- predict(g1, newx = X_test[, setdiff(colnames(X_test), monov)], s=maxL, type = 'response' )
  
  if(glmnetFam == "binomial" & dim(pred.response)[2] == 1){
    lbs <- levels(Y_train) # we have verified before that Y_train has to be a factor
#    lb.exist <- colnames(pred.response)
#    lb.missing <- setdiff(lbs, lb.exist)
# it seems pred.response may just have column name of '1' and igore the levels of Y_train    
    pred.response <- cbind(pred.response, (1 - pred.response))
    colnames(pred.response) <- rev(lbs) # c(lb.exist, lb.missing)
  }
  
  ## for survival, change the column back so that it is consistent across different algo.
  if(glmnetFam == "cox") {
    colnames(Y_test) <- c('col_surv','col_event') 
    }
  
  if(glmnetFam %in% c("multinomial", "binomial") ){
    pred.class <- predict(g1, newx = X_test[, setdiff(colnames(X_test), monov)], s=maxL, type = 'class')
    pred.out <- data.frame(Y_test,
                           pred.class,
                           pred.response,
                           stringsAsFactors = F)
    colnames(pred.out) <- c(colnames(Y_test),
                            "pred", paste0("predprob_", colnames(pred.response)))
  } else {
    pred.out <- data.frame(Y_test,
                           pred.response,
                           stringsAsFactors = F)
    colnames(pred.out) <- c(colnames(Y_test), "pred") ## assuming cox output only contains one column, and also name it 'pred'
  }

  if(is.list(g1$beta)) {
    b1 <- do.call(cbind, g1$beta)
  }else{
    b1 <- g1$beta
  }
  
  ## ToDo: temporarily use abs(x) > 0 to extract important features.
  ## need to use coefficient in the future. 
  b1 <- as.matrix(b1)
  features <- b1[apply(b1, 1, function(x) any(abs(x) > 0)),,drop=F]
  
  if(nrow(features) > 0) {
    features <- data.frame(row.names(features), features, row.names = NULL, stringsAsFactors = F)
  }else{
    features <- data.frame(matrix(rep(NA, ncol(b1) + 1), nrow=1))
  }

  colnames(features) <- c('feature', colnames(b1))  
  
  ## if all the coefficients are 0 (no important features)
#  if(nrow(features)==0){ features <- NA}
  
  list(features = features, 
       lambda = maxL, 
       pred = pred.out) 
}


rfeSRCCv3 <- function(X_train, Y_train, X_test, Y_test, sizes, seed, 
                      RFE_criteria = "permute", outputPrediction = 'chf'){
  require(prodlim)
  require(survival)
  
  set.seed(seed)
  
  fl <- names(X_train)
  df_vimp <- data.frame()
  df_pred <- data.frame()
  
  for(s in sizes){
    currfl <- tail(fl, s) 
    
    r1 <- randomForestSRC::rfsrc(Surv(col_surv, col_event) ~ ., 
                                 data = cbind(X_train[, currfl], Y_train), 
                                 importance = RFE_criteria, ## runSpec$RFE_criteria, currently fix it to permute
                                 seed = seed) 
    
    importance <- r1$importance
    stopifnot(all(sort(names(importance)) == sort(currfl)))
    importance <- sort(importance)
    fl <-  names(importance) ## sorted feature list
    
    ## vimp
    df_vimp_tmp <- data.frame(vimp = importance , 
                              feature = names(importance), 
                              size = s, 
                              stringsAsFactors = F, 
                              row.names = NULL)
    ## prediction
    pred <- predict(r1, newdata = X_test[, currfl, drop = F])
    
    ## survival prob
    pred2 <- NULL
    if ( nchar(outputPrediction) > 0 && outputPrediction %in% names(pred) ) {
      pred2 <- pred[[outputPrediction]]
      if (outputPrediction %in% 'chf'){pred2 <- exp( -pred2 )} # convert cumulative hazard to survival
      pred2 <- data.frame(pred2)
      rownames(pred2) <- rownames(X_test)
      pred2.times <- as.character( pred$time.interest )
      pred2 <- apply(pred2 , 1, paste, collapse = "," ) ## squeeze all columns to one single column
    }
    
    ## get brier score
    if( dim(X_test)[1] > 1 ) {
      bs1 <- pec::pec(object = list('randomforestSRC' = r1), 
                      formula = Surv(col_surv, col_event) ~ 1,  
                      data = cbind(X_test, Y_test), 
                      exact = TRUE, cens.model="marginal", splitMethod="none", B=0, verbose=F)
      ibrier <- pec::crps(bs1)
    } 
    
    df_pred_tmp <- data.frame(pred = pred$predicted, 
                              size = s,
                              BrierScore = ifelse( dim(X_test)[1] > 1, ibrier['randomforestSRC', 1], NA),
                              surv_prob = pred2,
                              stringsAsFactors = F, 
                              row.names=NULL)
    
    df_pred_tmp <- cbind(Y_test, df_pred_tmp)
    ## save surv_prob times, because different CV may have different predicted prob times
    ## will unlist this information during harvesting
    ## this is to make the workflow and output format consistent with xgboost
    df_pred_tmp$surv_prob.times = list(pred2.times) 
    
    df_vimp <- rbind(df_vimp, df_vimp_tmp)
    df_pred <- rbind(df_pred, df_pred_tmp)
  }
  
  return(list(df_vimp = df_vimp,
              df_pred = df_pred))
}