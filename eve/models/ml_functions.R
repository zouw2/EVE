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

extractBeta <- function(gObj, lambda){
  stopifnot('glmnet' %in% class(gObj))
  b0 <-   coef(gObj, lambda) 

  if(is.list(b0)) {
    b1 <- do.call(cbind, b0)
    colnames(b1) <- paste('vimp_', names(b0), sep='')
  }else{
    b1 <- b0
    stopifnot(ncol(b1) == 1)
    colnames(b1) <- 'vimp' #to work with reporting program
  }

  b1 <- as.matrix(b1)
  b1 <- b1[setdiff(row.names(b1), "(Intercept)"), ,drop=F]
  b1[apply(b1, 1, function(x) any(abs(x) > 0)),,drop=F]
}


glmnetCVwrapper2 <- function(X_train , Y_train, X_test, Y_test,
                             seed = 27519, 
                             glmnetFam="binomial", a1=1, 
                             nCv4lambda=10, lambdaSum=mean, 
                             runPairs=c(),
                             lambdaChoice = 'lambda.min',  w=1, ... ){
  
  set.seed(seed)
  stopifnot(is.matrix(X_train)); stopifnot(is.matrix(X_test))
  stopifnot(rownames(X_test) == rownames(Y_test))
  
  if (glmnetFam %in% c("multinomial", "binomial") ) {
    if (class(Y_train) %in% c('matrix','data.frame')) Y_train <- Y_train[, 1]
    #   if (class(Y_test) %in% c('matrix','data.frame'))  Y_test <-  Y_test[, 1]
    stopifnot( is.factor(Y_train) );
    #    stopifnot( is.factor(Y_test) );
    if (nlevels(Y_train) == 2) stopifnot(glmnetFam %in% 'binomial')
    if (nlevels(Y_train) > 2) stopifnot(glmnetFam %in% 'multinomial')
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
      w <- w[!idx2drop] ## drop the weights associated with those samples as well
    }
  }
  
  if (is.null(w)){ ## if no weight vector is given
    w = rep(1, dim(X_train)[1])
  }else{
    if(glmnetFam %in% c("multinomial", "binomial")){
    print('weight distribution by outcome levels')
    print(table(w,Y_train))
    }
  }
  
  print(paste('alpha=', a1))
  
  print('dim of training x')
  print(dim(X_train))
  
  f <- setdiff(colnames(X_train), findConstantVar(X_train))
  
  x1 <- X_train[, f]

  maxL <- lambdaSum( sapply(1:nCv4lambda, function(i) {
    set.seed(seed +100 + i)
    
    r1 <- cv.glmnet(x=x1, y=(function(yobj){
      if(glmnetFam =='cox') return(Surv(yobj)); 
      return(yobj) })(Y_train), family=glmnetFam, alpha = a1 , standardize =T, weights=w, ...)
    
    r1[[lambdaChoice]]
    
  }) , na.rm = T)

  
  g1 <- glmnet(x=x1, y = (function(yobj){
    if(glmnetFam =='cox') return(Surv(yobj)); 
    return(yobj) })(Y_train),  family = glmnetFam, alpha = a1,  standardize =T, weights=w, ...)

# features from g1
  features <- extractBeta(g1, lambda=maxL) 
  featureCol <- colnames(features)
 
  
# feature select with runPairs
  subsetFeatures <- c()
  if(length(runPairs) > 0){
    stopifnot(glmnetFam %in% "multinomial")
    
    pairwiseFeatures <-  lapply(runPairs, function(g){
      g <- unique(g)
      
      print(paste('select features between', paste(g, collapse = ',')))
      
      stopifnot(length(g) > 1)
      
      if(!all(g %in% levels(Y_train))){
        print( paste('not all levels from', paste(g, collapse = ','),'are found in this training data; no features will be selected from this comparison'))
        return(c())
      }
      
      sel2 <- as.character(Y_train) %in% g
      
      x2 <- x1[sel2, ]; x2[, setdiff(colnames(x2),  findConstantVar(x2))]
      y2 <- factor(as.character(Y_train[sel2])); 
      stopifnot(all(as.character(y2) %in% g));
      
      fam2 <- ifelse(length(g) ==2, "binomial", "multinomial")
      
      maxL2 <- lambdaSum( sapply(1:nCv4lambda, function(i) {
        set.seed(seed + 200 + i)
        r1 <- cv.glmnet(x=x2, y=y2, family=fam2, alpha = a1 , standardize =T, weights = w[sel2], ...)
        r1[[lambdaChoice]]
        
      }) , na.rm = T)
      
      g2 <- glmnet(x=x2, y = y2,  family = fam2, alpha = a1, standardize =T, weights = w[sel2], ...)
      
      row.names(extractBeta(g2, lambda=maxL2 ))
      
    })
    
    subsetFeatures <- unique(do.call(c, pairwiseFeatures))
  }
  
  subsetFeatures <- setdiff(subsetFeatures, row.names(features))
  
  
# final prediction  
  if(length(subsetFeatures) == 0){ # feature selection and outcome prediction will be the same model
    
  }else{ # refit a ridge regression and use that model to prediction. Note the coefficient values from this approach are not comparable to the standard glmnet flow
    print(paste('analyses in subgroups using a subset of labels adds', length(subsetFeatures), 'to plain glmnet selection of', nrow(features),'features'))
    f <- unique(c( row.names(features), subsetFeatures))
    
    maxL  <- lambdaSum( sapply(1:nCv4lambda, function(i) {
      set.seed(seed + 300 + i)
      r1 <- cv.glmnet(x=x1[, f,drop=F], y=Y_train, family= glmnetFam, alpha = 0  , standardize =T, weights=w, ...)
      r1[[lambdaChoice]]
      
    }) , na.rm = T)
    
    g1 <- glmnet(x=x1[, f,drop=F], y = Y_train,  family = glmnetFam, alpha = 0, standardize =T, weights=w, ...)
    
    features <- extractBeta(g1, lambda=maxL)
    
  }

  pred.response <- predict(g1, newx = X_test[, f,drop=F], s=maxL, type = 'response' )
    
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
    pred.class <- predict(g1, newx = X_test[, f,drop=F], s=maxL, type = 'class')
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


  # prepare feature matrix
  
  if(nrow(features) > 0) {
    features <- data.frame(row.names(features), features, row.names = NULL, stringsAsFactors = F)
  }else{
    features <- data.frame(matrix(rep(NA, ncol(features) + 1), nrow=1))
  }

  colnames(features) <- c('feature', featureCol)  
  
  ## if all the coefficients are 0 (no important features)
#  if(nrow(features)==0){ features <- NA}
  
  list(features = features, 
       lambda = maxL, 
       pred = pred.out) 
}


rfeSRCCv3 <- function(X_train, Y_train, X_test, Y_test, sizes, seed, 
                      RFE_criteria = "permute", outputPrediction = 'chf', ...){
  require(prodlim)
  require(survival)
  
  set.seed(seed)
  
  fl <- names(X_train)
  df_vimp <- data.frame()
  df_pred <- data.frame()
  
  for(s in sizes){
    currfl <- tail(fl, s) 
    
    r1 <- randomForestSRC::rfsrc(Surv(col_surv, col_event) ~ ., 
                                 data = cbind(X_train[, currfl,drop=F], Y_train), 
                                 importance = RFE_criteria, ## runSpec$RFE_criteria, currently fix it to permute
                                 seed = seed, ...) ## use mc.cores = 1 when running rfsrc in rstudio in rescomp
    
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
    pred2.times <- pred2 <- NA
    if ( nchar(outputPrediction) > 0 && outputPrediction %in% names(pred) ) {
      pred2 <- pred[[outputPrediction]]
      if (outputPrediction %in% 'chf'){pred2 <- exp( -pred2 )} # convert cumulative hazard to survival
      pred2 <- data.frame(pred2)
      rownames(pred2) <- rownames(X_test)
      pred2.times <- as.character( pred$time.interest )
      colnames(pred2) <- pred2.times
      ## break pred2 rows into list, each item is one row, 
      ## this is to help squeeze the probability prediction values into one single cell
      pred2 <- lapply(c(1:dim(pred2)[1]), function(x) pred2[x,]) 
      #pred2 <- apply(pred2 , 1, paste, collapse = "," ) ## squeeze all columns to one single column
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
                              #surv_prob = pred2,
                              stringsAsFactors = F, 
                              row.names=NULL)
    
    df_pred_tmp <- cbind(Y_test, df_pred_tmp)
    ## save surv_prob times, because different CV may have different predicted prob times
    ## will unlist this information during harvesting
    ## this is to make the workflow and output format consistent with xgboost
    df_pred_tmp$surv_prob <- pred2
    #df_pred_tmp$surv_prob.times = list(pred2.times) 
    
    df_vimp <- rbind(df_vimp, df_vimp_tmp)
    df_pred <- rbind(df_pred, df_pred_tmp)
  }
  
  return(list(df_vimp = df_vimp,
              df_pred = df_pred))
}


alignProb <- function(pred, timeNeeded){
  stopifnot(is.data.frame(pred$surv_prob[[1]]))
  timeAval <- colnames( pred$surv_prob[[1]] )
  stopifnot( all(sapply(pred$surv_prob, function(d) all(colnames(d) == timeAval))) ) # all elements have the same columns (probabilities from the same set of time points), as these are from the same cv
  
  p2 <- do.call(rbind, pred$surv_prob)
  stopifnot(all(row.names(p2) == row.names(pred)))
  
  colnames(p2) <- make.names(colnames(p2))
  p3 <- fillInMissingProbWithNeighbor(p2, cols=timeNeeded, method='LinearInterpolation')
  stopifnot(all(row.names(p3) == row.names(p2)))
  
  colnames(p3) <- gsub('^X','', colnames(p3))
  stopifnot(is.data.frame(p3))
  
  pred$surv_prob <- lapply(1:nrow(p3), function(x) p3[x,,drop=F]) 
  pred
}

fillInMissingProbWithNeighbor <- function(x, cols, pattern='^X[\\d\\.]+$', method = c('previousEstimate', 'LinearInterpolation' )){
  
  stopifnot(all(grepl('^X', colnames(x), perl=T)))
  cols <- paste0("X", cols)
  
  extra <- setdiff(colnames(x), cols)
  if(length(extra) > 0) print(paste('survival probs from the following timepoints are ignored in pre-validation brier score calculation', paste(extra, collapse = ',')))
  
  if(all(cols %in% colnames(x))) return(x[, cols, drop=F])
  
  stopifnot(! 'one' %in% colnames(x) )

  tofill <- setdiff(cols, colnames(x))

  if(length(tofill)> 0 && method == 'LinearInterpolation'){
    ta1 <- as.numeric( gsub('^X','', colnames(x), perl=T) )
    stopifnot(!any(is.na(ta1)))
    stopifnot(all(ta1 >= 0))
    
    fillin <- sapply(tofill, function(v, timeAvailable = sort( ta1 )){
      v <- as.numeric( gsub('^X','', v, perl=T) )

      stopifnot(!is.na(v))
      stopifnot(v >= 0)
      stopifnot(! v %in% timeAvailable )
      
      if(all(timeAvailable < v)){
        x.last <- x[, paste('X', tail(timeAvailable, 1), sep='')]
        #    colnames(x.last) <- paste0("X", v)
        return( unname(unlist(x.last)) ) # using the last time point
      }
      
      rightTime <- timeAvailable [  min( which(timeAvailable > v))    ]
      rightProb <- x[, paste('X', rightTime, sep='')]
      rightProb <- as.numeric(unlist(rightProb))
      
      if(all(timeAvailable > v)) {
        stopifnot(min(timeAvailable) > 0)
        leftProb <- rep(1, nrow(x))
        leftTime <- 0
      }else{
        leftTime <- timeAvailable [  max( which(timeAvailable < v) ) ]
        leftProb <- x[, paste('X', leftTime, sep='')]
        leftProb <- as.numeric(unlist(leftProb))
      }
      
      (leftProb * (rightTime -v)  + rightProb * ( v - leftTime ) )/(rightTime - leftTime)
      
      #this is way too slow
      #return( apply(cbind(leftProb, rightProb), 1, function(p){approx(x=c(leftTime, rightTime), y=p, xout=v )$y}) )
      
      
    })
    
    if(F){
      ## ToDo: temparaily fix
      ## There is a bug, for example, 
      ## if last time point is imputed, the column name will become X639.X639
      ## the temp fix is to remove any column names that has second ".X", then remove everything after
      fillin <- data.frame(t(unlist(fillin)))
      colnames(fillin) <- gsub("\\.X.*", "", colnames(fillin))
      #names(fillin) <- gsub("\\.X.*", "", names(fillin))
      #  colnames(fillin ) <- tofill
    }
    
    if(is.vector(fillin)) fillin <- matrix(fillin, nrow=1, dimnames=list(c(), names(fillin)))
    
    stopifnot(all(colnames(fillin) == tofill))
    x <- cbind(x, fillin)
    
  }
  
  if(method == 'previousEstimate'){
    stop('this option may be outdated for the current intended use of this function')
    fillin<- sapply(tofill, function(v, timeAvailable = as.numeric( gsub('^X','', setdiff(cols, tofill), perl=T) ) ){
      v <- as.numeric( gsub('^X','', v, perl=T) )
      stopifnot(!is.na(v))
      stopifnot(!any(is.na(timeAvailable)))
      
      stopifnot(! v %in% timeAvailable )
      timeAvailable <- sort(timeAvailable)
      
      if(all(timeAvailable > v)) return('one') # earlier than any observed event time, the survival prob will be set to 1
      i <- max( which(timeAvailable < v)  )
      v1 <- paste('X', timeAvailable[i], sep='')
      #browser()
      #    print(paste('fill in time', v, 'survival probability with those from time', v1))
      v1
    })
    
    if(! all( fillin %in% c('one', colnames(x)) )) {
      stop(paste('column', paste(setdiff(fillin, colnames(x)), collapse = ','), 'could not be found'))
    }
    
    
    x1 <- cbind(x, one=1)[, fillin, drop=F]
    colnames(x1) <- tofill # so even a column of survival probabilities of 1 is filled, the column name does not indicate this column is for time 0
    
    x <- cbind(x, x1)
  }
  
  
  
  stopifnot( all(cols %in% colnames(x)) )
  x.interpolate <- x[, cols, drop=F]

  return(x.interpolate)
}
