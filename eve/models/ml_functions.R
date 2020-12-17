vbind <- function(lst) {
  
  stopifnot(is.list(lst))
  n1 <- sapply(lst, is.null)
  if(any(n1)){
    lst <- lst[!n1]
    print(paste('removing', sum(n1),'elements with NULL values'))
    if(length(lst) == 0) return (NULL)
  }
  stopifnot(all(sapply(lst, is.vector)))
  stopifnot(all(sapply(lst, function(x) length(x) == length(names(x))))) 
  
  nam1 <- Reduce(union, lapply(lst, names))
  
  r1 <- matrix(NA, nrow=length(nam1), ncol=0, dimnames = list(nam1, c()))
  
  t( Reduce(cbind, lapply(lst, function(x) x[nam1]), init = r1) )
  
}

mean1 <- function(a){
  stopifnot(is.vector(a))
    a <- a[!is.na(a)]
    sum(a) / length(a)
}

harm.mean <- function(a){
    stopifnot(all(a > 0))
    1/mean1(1/a)
}

#both preds and labels use 0, 1, ... to indicate different levels. So the lowest level is 0, and additional levels add up.

# in the returning matrix, row status correponds to prediction, and column status for truth( labels)
create.xtab.nclass <- function(preds, labels, nclass=2 ) {
#modified from https://gist.github.com/khotilov/a0bd2f617e11811ede89287e42796ba5  
#  stopifnot( nclass == data.table::uniqueN(labels)) # the label vector must have all levels
  
  m1 <- matrix(tabulate(1 + preds + nclass*labels, nclass^2), ncol=nclass)
  stopifnot(nrow(m1) == nclass)
  stopifnot(ncol(m1) == nclass)
  m1
}

pr_scores <- function(cm, substitute0 = 0.5){
  # from https://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html#perclass
  # return value: element 1: precision of level 0, 1,2... ; element 2: recall of level 0,1,2
  diag = diag(cm)
  diag <- ifelse(diag == 0, substitute0, diag)
  sum_pred = apply(cm, 1, sum) 
  sum_label= apply(cm, 2, sum)
  
  return( list(precision = ifelse(sum_pred==0, substitute0/sum(sum_label), diag / sum_pred),   
               recall =   diag / sum_label))
}


# using harmonic means seem to have too much contribution by the rarer subgroup. evalution metrics never improve so early stopping occurs very early
#preds must be numbers
f1_harmonic2 <- function(preds, dtrain) {
  Actual = getinfo(dtrain, "label")
  lv_a <-  sort(unique(Actual))
  nc <- length(lv_a)
  stopifnot(nc == 2)
  if(! all(as.vector(preds) %in% lv_a) ) preds <- as.integer(preds >= 0.5)
  stopifnot(all( as.vector(preds)  %in% lv_a))
  value <- harm.mean(unlist(pr_scores(create.xtab.nclass(preds, Actual,nclass=nc))))
  return(list(metric = "f1_harmonic2", value = value))
}

# log scale mean is equivalent to geometric mean. it does not help with orr prediction in oak either.
#preds must be numbers
f1_meanlog2 <- function(preds, dtrain) {
  Actual = getinfo(dtrain, "label")
  lv_a <-  sort(unique(Actual))
  nc <- length(lv_a)
  stopifnot(nc == 2)
  if(! all(as.vector(preds) %in% lv_a) ) preds <- as.integer(preds >= 0.5)
  stopifnot(all( as.vector(preds)  %in% lv_a))
  value <- mean1(log(unlist(pr_scores(create.xtab.nclass(preds, Actual,nclass=nc)))))
  return(list(metric = "f1_meanlog2", value = value))
}


#preds must be numbers
f1_arithmetic2 <- function(preds, dtrain) {
  Actual = getinfo(dtrain, "label")
  lv_a <-  sort(unique(Actual))
  nc <- length(lv_a)
  stopifnot(nc == 2)
  if(! all(as.vector(preds) %in% lv_a) ) preds <- as.integer(preds >= 0.5)
  stopifnot(all( as.vector(preds)  %in% lv_a))
  value <- mean1( unlist(pr_scores(create.xtab.nclass(preds, Actual,nclass=nc)))) 
  return(list(metric = "f1_arithmetic2", value = value))
}

f1_1side_2 <- function(preds, dtrain) {
  Actual = getinfo(dtrain, "label")
  lv_a <-  sort(unique(Actual))
  nc <- length(lv_a)
  stopifnot(nc == 2)
  if(! all(as.vector(preds) %in% lv_a) ) preds <- as.integer(preds >= 0.5)
  stopifnot(all( as.vector(preds)  %in% lv_a))
  pr <- pr_scores(create.xtab.nclass(preds, Actual,nclass=nc))
  value <- min( sapply( 1:nc, function(i) harm.mean(c( pr[['precision']][i], pr[['recall']][i])) ) ) 
  return(list(metric = "f1_1side_2", value = value))
}

print1 <- function(x, ...) print(paste(names(x),':', x, collapse = ' '), ...)

xgbCVwrapper <- function(X_train, Y_train, X_test, Y_test , ft_seqs,  weight_train, params, spec, min_n_estimators=5 ){
  
  if(!is.null(spec$nthread) && spec$nthread > 1 && !('nthread' %in% names(params))){
    params[['nthread']] <- spec$nthread 
    print(paste('insert nthread of',spec$nthread ,'to params'))
  }else{
    params[['nthread']] <- 8
    print('using default nthread of 8')
  }
  
  internal_eval_metric <- c('error','rmse','logloss','mlogloss','auc','aucpr','merror', 'ndcg')
  max_not_min <- grepl('auc|f1', params$eval_metric)
  
  per_cv_seed <- (spec$'localSeed' - spec$seed_base) * (nrow(X_train) + nrow(X_test) + 5) + spec$curr_cv #this seed should be different per cv per split/seed
  
  stopifnot(all(colnames(X_train) == colnames(X_test)))
  
  # allow users to replace default params specified in XGBoost.r, only if a single value is specified
  if(length(spec$tune_params) > 0 ) {
    len <- sapply(spec$tune_params, length)
    for (v in names(len)[len==1]){
      params[[v]] <- spec$tune_params[[v]]
      print(paste('assign', v,'as', spec$tune_params[[v]]))
    }
  }

  if(max(ft_seqs) < ncol(X_train)){
    print('add a full feature step to the rfe')
    ft_seqs <- c(ncol(X_train), ft_seqs)
  }
  
  rfe_grid <- ft_seqs
  
  if( !(is.null(spec$pct_feature_for_grid_search) || is.na(spec$pct_feature_for_grid_search)) && 
    is.numeric(spec$pct_feature_for_grid_search) && 
    length(ft_seqs) > length(spec$pct_feature_for_grid_search)){
    r <- ncol(X_train) * spec$pct_feature_for_grid_search
    rfe_grid <- sapply(r, function(v) ft_seqs[which.min(abs(v - ft_seqs))])
  }

  
  
  library(xgboost)
  # from https://pages.github.roche.com/RWDScodeshare/RAADC_2018site/tutorials_xgboost/

    
  stopifnot(is.vector(Y_train))
  stopifnot(is.vector(Y_test))

  # Create folds within X_train. folds are the same across parameter grid (following roche tutorial)
  # so far I put it above RFE
  set.seed(spec$'localSeed')
  folds1 <- stratFold ( outcome =  Y_train , fam = spec$family, num_CV = spec$nCv4lambda)

  features <- colnames(X_train)

  vimpL <- gridL <- predL <- list()
  params2 <- params
  
  for (k in ft_seqs){
  
    # update the following object during RFE
    top_fts = head( features, k) # features vector will be modified during RFE, top_fts holds the current features used in training
  
    stopifnot(all(top_fts %in% colnames(X_train)))
    dtrain <- xgb.DMatrix(
      data =  X_train[, top_fts] , 
      label = Y_train,  # not do as.matrix as suggested originally
      weight = weight_train
    )
  
  
  if(length(spec$tune_params) > 0 && k %in% rfe_grid){
    
    stopifnot(! 'eval_metric' %in% names(spec$tune_params)) # Not sure how to tune among different evaluation metric as they have different scale
    e1 <- expand.grid(spec$tune_params, KEEP.OUT.ATTRS =F, stringsAsFactors=F)
    paraList <- split(e1, seq(nrow(e1)))
    paraList <- lapply(paraList, function(x) as.list(x))
    paraList <- lapply(paraList, function(x) c(x, params[setdiff(names(params), names(x))]))
    
    if((!is.null(spec$max_num_grid)) && (!is.na(spec$max_num_grid)) && spec$max_num_grid > 0){ # if a user provides max_num_grid
    set.seed( per_cv_seed ) # so if there are many combinations, different CVs within the same seed will try different combinations
    paraList <- sample( paraList ) # reshuffle it
    }else{ # if a user does not provide max_num_grid
      spec$max_num_grid <- length(paraList)
    }
#============grid search
    print(paste('using', ifelse(max_not_min, 'max', 'min'),  params$eval_metric, 'to select the best in', min(length(paraList), spec$max_num_grid),'parameter combinations out of', length(paraList),'all possibilities at feature step', length(top_fts))) 
    
    start.time <- Sys.time()
    
    grid_best <- do.call(rbind, lapply( 1: min(length(paraList), spec$max_num_grid), function(i) {
      print('parameters being evaluated:')
      print1(paraList[[i]])
      if(params[['eval_metric']] %in% internal_eval_metric ){
        cv1 <- xgboost::xgb.cv(
          params = paraList[[i]],
          data = dtrain,
          nrounds = ifelse(!(is.null(spec$n_estimators)||is.na(spec$n_estimators)), spec$n_estimators, 1000),
          folds = folds1, # so the cv split will be the same across the entire unit run
          prediction = F,
         # verbose = i ==1,
         verbose=F,
          print_every_n = 100,
          metrics = params[['eval_metric']],
          maximize =  max_not_min,
          
          early_stopping_rounds = ifelse(is.null(spec$early_stopping_rounds), 
                round(2/ifelse(is.null(paraList[[i]]$eta), params$"eta", paraList[[i]]$eta ) ) , spec$early_stopping_rounds)      )
        }else{ # use user defined feval
        cv1 <- xgboost::xgb.cv(
          params = paraList[[i]][setdiff(names(paraList[[i]]),'eval_metric')],
          data = dtrain,
          nrounds = ifelse(!(is.null(spec$n_estimators)||is.na(spec$n_estimators)), spec$n_estimators, 1000),
          folds = folds1, # so the cv split will be the same across the entire unit run
          prediction = F,
          # verbose = i ==1,
          verbose=F,
          print_every_n = 100,
          feval = match.fun(params[['eval_metric']]),
          maximize =  max_not_min,
          
          early_stopping_rounds = ifelse(is.null(spec$early_stopping_rounds), round(2/ifelse(is.null(paraList[[i]]$eta), params$"eta", paraList[[i]]$eta ) ), spec$early_stopping_rounds)  
        )
      }
      
      el <- as.data.frame(cv1$evaluation_log)
      
      cv1$score <- el[el$iter == cv1$best_iteration, grep('test.*mean$', colnames(el), perl=T, value=T)]
      
      cv2 <- data.frame(size = cv1$nfeatures, score = cv1$score, n_estimators = cv1$"best_ntreelimit", cv1$params[setdiff(names(cv1$params), c( 'weight', 'silent'))], stringsAsFactors = F)
      
#      cv2$weight <- list(unique(cv1$params$weight))
      cv2
    }) )
    
 
    gridL[[as.character(k)]] <- grid_best
    gridL[[as.character(k)]]$cv <- spec$curr_cv
    end.time <- Sys.time()
    
    print(paste('grid search takes: ', round(difftime(end.time , start.time, units ="hours"),2),"hours")) 
    
    # end of grid search
    
    params2 <- as.list( grid_best[which.max( grid_best$score * ifelse(max_not_min, 1, -1)  ),  ] )
 } 
  
  print( paste( 'parameter at step',k))   
  print1(unlist( params2[setdiff(names(params2), c('size','score','weight','silent'))] ))
  
  set.seed(spec$'localSeed')
   
  nd1 <- ifelse('n_estimators' %in% names(params2), params2$n_estimators, 5/params2$eta) #https://www.slideshare.net/OwenZhang2/tips-for-data-science-competitions slide 14
  
  if(nd1 < min_n_estimators) nd1 <- min_n_estimators
  
   xgb.model <- xgboost::xgb.train(
        params = params2[setdiff(names(params2), c('size','score','n_estimators','weight','silent','eval_metric'))],
        data = dtrain,
        nrounds = nd1
      )
 
     #get variable importance
   
   vimp <- as.data.frame ( xgb.importance(model = xgb.model) )
   colnames(vimp) <- tolower(colnames(vimp))
   
   stopifnot(spec$RFE_criteria %in% colnames(vimp))
   
   vimp <- vimp[order(vimp[, spec$RFE_criteria] * -1), ]
   # when the max_depth is low and the n_estimates is not too big, it is possible that some features are not used in the final model and vimp matrix does not contain the variable. Need figure out a univariate way to sort genes
   
   if(nrow(vimp) < length(top_fts)){
     set.seed( per_cv_seed )
     features <- c(vimp$feature, sample(setdiff(top_fts, vimp$feature))) # so here features vector is shorter than full feature list, but has the same length of top_fts
   }else{
     features <- vimp$feature
   }
   
   # add columns to vimp to confirm to the reporting standard
   
   colnames(vimp) <- gsub('frequency','weight', colnames(vimp))
   vimp$size <- length(top_fts)
   vimp$cv  <-  spec$curr_cv
   vimp$cover <- NULL
   
   vimpL[[as.character(k)]] <- vimp
   
   # Predict on the standalone test set
   stopifnot(all(top_fts %in% colnames(X_test)))
   
   if (spec$family %in%  "gaussian") {
     pred <-      stats::predict(xgb.model, X_test[,top_fts,drop=F])
     xgb_prediction <- data.frame(pred,  stringsAsFactors = F)
   }
   
   if( spec$family %in% 'binomial') {
     predprob_1 = stats::predict(xgb.model, X_test[,top_fts,drop=F]) # the numbers we get are probabilities that a datum will be classified as 1
     stopifnot(length(predprob_1) == nrow(X_test))
     predprob_0 <- 1- predprob_1
     xgb_prediction <- data.frame(pred = ifelse(predprob_1 > 0.5, 1,0), predprob_0, predprob_1,  stringsAsFactors = F)
     }   
   
   if(  spec$family %in% "multinomial" ){
     predprob <- stats::predict(xgb.model, X_test[,top_fts,drop=F], reshape=T)
     stopifnot(ncol(predprob) == params2$num_class)
     colnames(predprob) <- paste('predprob_', 0:(params2$num_class-1), sep='')
     xgb_prediction <- data.frame(predprob)
     xgb_prediction$pred <- apply(predprob, 1, which.max) -1
   }
   
   xgb_prediction[[spec$label_name]] <- Y_test
   xgb_prediction[['size']] <- length(top_fts)
   xgb_prediction[['cv']] <- spec$curr_cv
   xgb_prediction[[spec$sample_ID]] <- row.names(X_test)
      
   predL[[as.character(k)]] <- xgb_prediction
  }
  
  list(do.call(rbind, vimpL) , do.call(rbind, gridL) , do.call(rbind, predL))
} 


stratFold <- function( outcome, fam, num_CV){
  if (fam %in% 'cox') {
    if( all( c('time','status') == colnames(outcome) ) ) colnames(outcome) <- c('col_surv', 'col_event')
    
    stopifnot(all(c('col_surv', 'col_event') %in% colnames(outcome)))
    
    return(  caret::createFolds(outcome[, "col_surv"] * (as.integer(outcome[, "col_event"]) - 0.5), 
                                 k = num_CV, 
                                 list = TRUE, returnTrain = FALSE) )
  }
  
  stopifnot(is.vector(outcome) || is.factor(outcome))
  
  if (fam %in% c("multinomial", "binomial")) {
    cvList <- caret::createFolds(as.factor(outcome), 
                                 k = num_CV, 
                                 list = TRUE, returnTrain = FALSE)
  }
  
  if (fam %in% c("gaussian")) {
    stopifnot(is.numeric(outcome))
    cvList <- caret::createFolds(outcome, 
                                 k = num_CV, 
                                 list = TRUE, returnTrain = FALSE)
  }
  
  cvList
}

decideRFEseq <- function(RFE_step  , ft){
#  stopifnot( length(RFE_step) == 1 && RFE_step >= 0 )
  if( is.null(RFE_step) || is.na(RFE_step) || nchar(as.character(RFE_step))==0){
    print('no RFE specified, just use full feature set')
    return(ft)
  }

  if( !is.numeric(RFE_step) ){
      print(paste('could not interpret the input RFE_step:', paste(RFE_step, collapse=',')))
      print('just use full feature set')
      return(ft)
  }
  
## if user provides a vector of steps, just return as is
if( length(RFE_step) > 1 ) return( RFE_step )

## if user provides a single value
  if (RFE_step == 0) return(ft)

  if(is.integer(RFE_step)){
    return( seq(ft, RFE_step,  RFE_step*-1) ) ## use fixed step_size
  }else{
    stopifnot(RFE_step > 1)
    step_size <-  min_step <- floor(RFE_step) # get the integer part as the min step size
    RFE_step <- RFE_step - min_step # get the fraction part as the fract of features to discard
    sizes <- c()
        
    while (step_size >= min_step){
      step_size <- round(RFE_step*ft)
      sizes <- c(sizes, ft) 
      ft <- ft - step_size
    }
 
    if(tail(sizes,1) > min_step * 2) {
      sizes <- c(sizes, seq(from= tail(sizes,1) - min_step, to = min_step, by=min_step*-1))
    }
    return(sizes)
  }
}

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
 #   colnames(b1) <- paste('vimp_', names(b0), sep='')
     colnames(b1) <- paste('coef_', names(b0), sep='')
    
  }else{
    b1 <- b0
    stopifnot(ncol(b1) == 1)
#    colnames(b1) <- 'vimp' #to work with reporting program
    colnames(b1) <- 'coef'  #to work with reporting program
  }

  b1 <- as.matrix(b1)
  b1 <- b1[setdiff(row.names(b1), "(Intercept)"), ,drop=F]
  b1[apply(b1, 1, function(x) any(abs(x) > 0)),,drop=F]
}




tuneLassoParam <- function (n, lsum, nfolds, x , y , fam, alpha , weights , ... ){
  # create folder id
  fold2  <- lapply(1:n, function(j){
    f1 <- stratFold(y, fam, nfolds)
    stopifnot(!any(duplicated(unlist(f1))))
    
    r1 <- rep(NA, length(unlist(f1)))
    for (i in 1:length(f1)) r1[f1[[i]]] <- i
    stopifnot(!any(is.na(r1)))
    stopifnot(all(sort(unique(r1)) == 1:nfolds))
    r1
  })
  
  if(fam == 'cox') y <- Surv(y)
  
  res <- lapply(alpha, function(a){
    
    r2 <- t(sapply(1:n, function(i) { 
      
      cv1 <- cv.glmnet(x=x, y=y, family = fam, foldid = fold2[[i]], alpha=a,  standardize =T, weights=weights, ...)
      
      return(c(alpha= a, lambda.min=cv1$lambda.min, lambda.1se = cv1$lambda.1se, cvm.min=min(cv1$cvm)))
      } )) } )

  } 
  
   




tuneLassoParam_averageCVM <- function (n, lsum, nfolds, x , y , fam, alpha , weights , ... ){
  # this approach to average cvm and find the (alpha, lambda) pair that minimized the averaged cvm does not work because the lambda sequences seem to different sometime even if I provided the lambda seqeunces and seed
  
  fold2  <- lapply(1:n, function(j){
    f1 <- stratFold(y, fam, nfolds)
    stopifnot(!any(duplicated(unlist(f1))))
    
    r1 <- rep(NA, length(unlist(f1)))
    for (i in 1:length(f1)) r1[f1[[i]]] <- i
    stopifnot(!any(is.na(r1)))
    stopifnot(all(sort(unique(r1)) == 1:nfolds))
    r1
    })
  
  if(fam == 'cox') y <- Surv(y)
  
  res <- lapply(alpha, function(a){
    
    r1 <- cv.glmnet(x=x, y=y, family = fam, foldid = fold2[[1]], alpha=a, standardize =T, weights=weights, ...)
    lseq <- r1$lambda #the labmda sequence is different for different alpha
    
    r2 <- lapply(2:n, function(i) { 
      set.seed(29517)
      cv.glmnet(x=x, y=y, family = fam, foldid = fold2[[i]], alpha=a, lambda=lseq, standardize =T, weights=weights, ...) } )
    
    s1 <- sapply(1:length(r2), function(i){ stopifnot( all( lseq == r2[[i]]$lambda ) )})
    
    r2[[length(r2)+1]] <- r1
    
    r2cvm <- t(sapply(r2, function(x) x$cvm)) # expecting a matrix (more than 1 alpha) with rows for different alpha and columns for different lambda
    
    
    stopifnot(nrow(r2cvm) == n)
    stopifnot(ncol(r2cvm) == length(lseq))
    
    return(list(lambda = lseq, cvm = apply(r2cvm, 2, lsum)))
    
  })
  
  if(length(alpha) == 1){
    stopifnot(length(res) == 1)
    with(res[[1]], stopifnot(is.numeric(cvm)))
    with(res[[1]], stopifnot(is.numeric(lambda)))
    return(list(lambda=with(res[[1]], lambda[which.min(cvm)]), alpha=alpha))
  }
  
  lm <- do.call(rbind, lapply(res, function(x) x$lambda))
  cvm <- do.call(rbind, lapply(res, function(x) x$cvm))
  print('initial error metrics')
  print(cvm[, 1:4])
  
  stopifnot(all(dim(lm) == dim(cvm)))
  stopifnot(nrow(cvm) == length(alpha))
  stopifnot(ncol(cvm) == length(res[[1]]$lambda))
  best <- arrayInd(which.min(cvm), dim(cvm))
  return( list( lambda = lm[best[1], best[2]], alpha = alpha[best[1]]))
}

glmnetCVwrapper2 <- function(X_train , Y_train, X_test, Y_test,
                             seed = 27519, 
                             glmnetFam="binomial", a1=1, 
                             nCv4lambda=10, lambdaSum=mean, 
                             runPairs=c(),
                             lambdaChoice = 'lambda.min',  w=1, usePerCVlambda=F, ... ){
  
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
    idx2drop <- is.na( Y_train[, "time"] ) | Y_train[, "time"] <= 0
    if(any(idx2drop)){
      print(paste("Drop", sum(idx2drop), "samples with survival time < 0 or missing."))
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
  
  print(paste('alpha=', a1, collapse = ','))
  
  print('dim of training x')
  print(dim(X_train))
  
  f <- setdiff(colnames(X_train), findConstantVar(X_train))
  
  x1 <- X_train[, f]

  a2 <- a1
  importance <- NULL
  
  set.seed(seed + 101 )
  
  if(length(a2) == 1){ # the previous way of only tuning lambda. It was commented out earlier, but wei brought this back when incorporating nextdoor analysis
    tuneResults <- lapply(1:nCv4lambda, function(i) {
  #    set.seed(seed +100 + i)
      print(paste('nested cv', i))
      
      
      r1 <- cv.glmnet(x=x1, y=(function(yobj){
        if(glmnetFam =='cox') return(Surv(yobj)); 
        return(yobj) })(Y_train), keep=T, family=glmnetFam, alpha = a2 , standardize =T, weights=w, ...)
      
      vimp <- NULL
      if(usePerCVlambda && require('nextdoor') && length(runPairs) == 0 && ! glmnetFam %in% c("multinomial")){  
        n1 <-  nextdoor.glmnet(x=x1, y=(function(yobj){
          if(glmnetFam =='cox') return(Surv(yobj)); 
          return(yobj) })(Y_train), cv_glm =r1, nams= colnames(x1), family=glmnetFam,  glmnet_alpha = a2, standardize =T , selectionType = ifelse(lambdaChoice == "lambda.1se", 1, 0 ), pv=F, score=F, trace = F) 
        
        vimp = unlist(n1$worsen)
        stopifnot(length(vimp) == length(n1$worsen)) # assume every element of worsen is a scalar
    #    vimp =  matrix(vimp, nrow=1, dimnames = list(c(),names(vimp)) )
        
      }

      list(lambda = r1[[lambdaChoice]], vimp=vimp, cv =  (function(){
        if(usePerCVlambda) return(NA)
        r1 # return the CV object for later next door analysis
        })() )
      
    })
    
    maxL <- lambdaSum( sapply(tuneResults, function(x) x[['lambda']])  , na.rm = T)
    
    if((!usePerCVlambda) && require('nextdoor')){ # override the default behavior of nextdoor (using getIndex() to find features to evaluate
      tuneResults <- lapply(tuneResults, function(tr){
        
        n1 <-  nextdoor.glmnet( x=x1, y=(function(yobj){
          if(glmnetFam =='cox') return(Surv(yobj)); 
          return(yobj) })(Y_train), cv_glm =tr[['cv']], nams= colnames(x1), family=glmnetFam,  glmnet_alpha = a2, standardize =T , selectionType = ifelse(lambdaChoice == "lambda.1se", 1, 0 ), pv=F, score=F, trace = F, sumLambda = maxL) 
        
        vimp = unlist(n1$worsen)
        stopifnot( length(vimp) == length(n1$worsen) )

        list(lambda = tr[['lambda']], vimp=vimp)
        
        })
      
    }
    
    
    
    if(!is.null(tuneResults[[1]][['vimp']] ) ) {
      importance <- vbind( lapply(tuneResults, function(x) x[['vimp']] ) ) 
      
      if(!is.null(importance)) {
        importance[is.na(importance)] <-  min( apply(importance, 2, min, na.rm=T), 0)
        importance <- colMeans(importance)
      }
    }
  }else{ # allow tuning alpha
  
  #  set.seed(seed + 101 )
    tuneResults <- tuneLassoParam(n=nCv4lambda, lsum =lambdaSum, nfolds=10, x=x1, y= Y_train, fam=glmnetFam, alpha=a1, weights=w, ... )
    
    tuneSummary <- t( sapply(tuneResults, function(x) {
        stopifnot(length(unique(x[, 'alpha'])) == 1)
      return(c(alpha = unname(x[1, 'alpha']), lambda=lambdaSum(x[, lambdaChoice]), cvm.min= mean(x[, 'cvm.min'])))
      }))
    
    print('tuning summary')
    print(tuneSummary)
    
    sel <- which.min(tuneSummary[,'cvm.min'])
    maxL <- tuneSummary[sel,'lambda']
    a2 <-   tuneSummary[sel,'alpha']
  }
  
  
  g1 <- glmnet(x=x1, y = (function(yobj){
    if(glmnetFam =='cox') return(Surv(yobj)); 
    return(yobj) })(Y_train),  family = glmnetFam, alpha = a2,  standardize =T, weights=w, ...)

# features from g1
  features <- extractBeta(g1, lambda=maxL) 
#  featureCol <- colnames(features)
 
  
# feature select with runPairs
  subsetFeatures <- c()
  if(length(runPairs) > 0){
    stopifnot(glmnetFam %in% "multinomial")
    if(length(a1) > 1) stop('alpha tuning has not been implemented for runPairs function')
    
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
    if(!is.null(importance)) {
      stopifnot(length(intersect( features[, 1] , names(importance) )) > 0) # importance may be calculate based on a more stringent lambda than maxL (determined by nextdoor::getIndex. so the length of importance can be shorter than features based on maxL
      features$vimp <- importance[features[, 1]]
    }  
  }else{
    features <- data.frame(matrix( NA, ncol = ncol(features) + 1, nrow=0))
  }

#  colnames(features) <- c('feature', featureCol)  
  colnames(features)[1] <- 'feature'
  
  ## if all the coefficients are 0 (no important features)
#  if(nrow(features)==0){ features <- NA}
  
  list(features = features, 
       lambda = maxL, alpha = ifelse(length(a1)==1, a1, a2),
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
    
    missed_feature <- setdiff( currfl, names(importance))
    if(length(missed_feature) > 0) {
      print(paste('there are', length(missed_feature), 'included in the input X matrix but excluded from importance vecor, for example:', paste(head(missed_feature, 10), collapse = ',')))
    }
    
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
    ta1 <- as.numeric( gsub('^X','', gsub('e\\.(\\d+)', 'e-\\1',  colnames(x), perl=T), perl=T) )
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
