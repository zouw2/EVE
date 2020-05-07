round.signif <- function(x, p)
{
  ifelse(abs(x)>=1, round(x, p), signif(x, p))
}


if( ! (R.version$major >= 3 && R.version$minor >= 5.1 )) {
  stop(paste('reporting needs R version >= 3.5.1'))
}

parseMajorMinor <- function(x){
  r <- as.integer(strsplit(x, split='.', fixed=T)[[1]])
  stopifnot(length(r) == 3)
  stopifnot(!any(is.na(r)))
  
  return(c('major'=r[1], 'minor'=r[2],'patch'=r[3]))
  }

myScatter <- function(ds, x, y, noise=c(x=0,y=0), col.var='', corMethod='', pch.var='', nCol=4, a=1, lowCol='green', midCol='black', highCol='red', manualCol=c(),  main=''){
  stopifnot(is.data.frame(ds))
  
  stopifnot(all(c(x,y) %in% colnames(ds)))
  
  
  
  names(noise) <- c(x, y)
  
  for( v in names(noise)) {
    if(noise[v] > 0)  ds[, v] <- ds[, v] + rnorm(n=nrow(ds), mean=0, sd=noise[v])  
  }
  
  ds <- ds[!(is.na(ds[, x]) | is.na(ds[, y])),,drop=F]
  
  if( nrow(ds) == 0 ) stop(paste('no observations left with both non-missing', x, 'and', y))
  
  p1 <- ggplot(ds, aes_string(x=x, y=y) )
  
  if( col.var %in% colnames(ds) ) {
  
    if(is.factor(ds[, col.var])){
      col.f.var <- col.var
      if(length(manualCol) > 0) {
        col1 <- manualCol
      }else{
        col1 <- gplots::colorpanel(n=nlevels(ds[[col.f.var]]), low=lowCol, mid=midCol, high=highCol)
      }  
    }else{
      col.f.var <- paste(col.var,'F', sep='_')
      stopifnot(! col.f.var %in% colnames(ds) )
    brk1 <- quantile( ds[, col.var], probs= seq(0.01, 0.99, length.out=nCol), na.rm=T )
    brk1 <- unique(c(min(ds[, col.var]), brk1, max(ds[, col.var])))
    
    ds[[col.f.var]] <- cut(ds[, col.var], breaks=brk1, include.lowest=T)
    if(length(manualCol) > 0) {
      col1 <- manualCol
    }else{
      col1 <- gplots::colorpanel(n=length(brk1) - 1, low=lowCol, mid=midCol, high=highCol)
    }
    
    }
    
    names(col1) <- levels(ds[[col.f.var]])
    
    if( pch.var %in% colnames(ds) ) {
      p1 <- p1 + geom_point(data=ds, aes_string( col=col.f.var, pch=pch.var), alpha=a) + scale_colour_manual( values=col1, name = paste(col.var,'bins') )
    }else{
      p1 <- p1 + geom_point(data=ds, aes_string( col=col.f.var ), pch=16,     alpha=a) + scale_colour_manual( values=col1 )
    }
  }else{
    if( pch.var %in% colnames(ds) ) {
      p1 <- p1 + geom_point(data=ds, aes_string(pch=pch.var), alpha=a) 
    }else{
      p1 <- p1 + geom_point( alpha=a)  
    }
    
     
  }
  
  main <- paste(main, 'n=', nrow(ds))
  
  if( corMethod %in% c("pearson", "kendall", "spearman")) {
    main <- paste(main, corMethod,'=', round.signif( cor(ds[, x], ds[, y], method=corMethod), 2))
  }
  
  p1 + ggtitle(main)

  # if( pch.var %in% colnames(ds) ) {
  #   ggplot(ds, aes_string(x=x, y=y)) + geom_point(aes_string(col=col.f.var, pch=pch.var)) + scale_colour_manual( values=col1, name = paste(col.var,'bins') )
  # }else{
  #   ggplot(ds, aes_string(x=x, y=y)) + geom_point(aes_string(col=col.f.var), pch=16) + scale_colour_manual( values=col1 )
  # }
  
  
}

countUniqueFeatures <- function(vimpIn, size ) {
  stopifnot(all( c('cv','seed') %in% colnames(vimpIn)))
  
  if ('size' %in% colnames(vimpIn)){
    if( missing(size) || is.null(size)|| is.na(size)) size <- max(vimpIn$size)
  
    vimpIn <-  vimpIn[vimpIn[, 'size'] == size, ,drop=F]
  }else{
    size <- 'current'
  }

  print(paste('there are', length(unique(vimpIn$feature)), 'unique features used from the', size,'feature set'))

  with(vimpIn, tapply(feature, list(seed,cv), function(x) {
 
    stopifnot( !any(duplicated(x)) )
    length(x)
    }))
}

extractParameterFromFileName <- function(fileName, pos, sep1='_'){
  fileName <- gsub('\\.[^.]*$', '', fileName, perl=T) #remove one file extension
  gsub('\\D','', strsplit(fileName, split=sep1)[[1]][pos], perl=T)
}

# this function checks for any missing cv splits
# in the file name, seed and cv numers are separated by _
checkOutput <- function(path, pattern = paste( '\\.rdata', sep=''), posSeed=0, posSplit = 0 ){
  # assumptions: no folder structure under the path; 
  newfiles <-  list.files(path=path)
  newfiles <- newfiles[grep(pattern, newfiles, perl=T)]
  print(paste('finding', length(newfiles),'input files  ' ))

  f1 <- newfiles[1]
  if(posSeed > 0)    print(paste('from filename', f1,'the seed value is extracted as', as.integer( extractParameterFromFileName(f1, posSeed) )))
  if(posSplit > 0)  print(paste('from filename', f1, 'the CV split number is extracted as', as.integer( extractParameterFromFileName(f1, posSplit) )))
  
  r2 <- t( sapply(newfiles, function(f){
    seed <- split <- NA
  if(posSeed > 0)    seed=as.integer( extractParameterFromFileName(f, posSeed) )
  if(posSplit > 0)  split=as.integer( extractParameterFromFileName(f, posSplit) )
    c(seed, split)
  }) )
  
  
  colnames(r2) <- c('seed', 'cv')
  maxCV <- max(r2[, 'cv'])
  seed1 <- sort(unique(r2[,'seed']))
  print( paste('I think this set of CV has', length(unique(r2[, 'seed'])),'seeds ranging from', seed1[1],'from', tail(seed1, 1),'with', maxCV,'folder CV') )  
  
  if(!all(diff(seed1) == 1)){ #seeds discontinues
    print('discontinuous seed regions:')
    for(i in which(diff(seed1) != 1)) print(seed1[i:(i+1)]) 
  }
  
  cv_count <- tapply(r2[,'cv'], r2[,'seed'], length)
  
  if(any(cv_count != maxCV)){
    print(paste('the following seeds have fewer than', maxCV,'CV'))
    print(cv_count[cv_count!= maxCV])
  }
  
}  
  
#' Read a single object from .rdata 
#' R based ML engines save outputs as .rdata, therefore,
#' we need to be able to parse those output
#' 
#' @param filepath path to .rdata
#' @param objName either df_pred or df_vimp
#' 
read_rdata <- function(filepath, objName){
  new_env <- new.env(parent = emptyenv())
  load(filepath, envir=new_env)
  
  if(exists(objName, envir=new_env)) {
    df <- new_env[[objName]]
  }else{
    warning(paste0(objName, " does not exist"))
  }
}

#' Parsing ML output files. 
#' This version is for harvesting one seed's CVs from different jobs.
#' The 'gsub' function is different from harvesting one seed's CVs in a single job.
#' 
#' @param project_home path to the project_home foler
#' @param project_name name of the project folder
#' @param type indicate the parser to parse 'error', 'vimp', 'gridsearch', or 'rdata' which contains df_pred and df_vimp objeects
#' @param objName (optional) For read_rdata; and possible values are df_pred or df_vimp
#' @param requireCPrevalidation by default, this program tries to exclude seeds with missing cross validation runs to enable prevalidation
#' 
getResults <- function(project_home, project_name, type, objName = NULL, requireCPrevalidation=T){
  Path2Results <- paste0(project_home, "/results/", project_name)
  #print(paste("reading files in", Path2Results))
  
  newfiles <- list.files(path=Path2Results)
  idx <- grepl(type, newfiles)
  idx_datatype <- grepl("\\.rdata$|\\.csv$", newfiles) ## only collect either .rdata or .csv files in that folder
  
  files <- newfiles[idx & idx_datatype] ## only collect either .rdata or .csv files in that folder
  #files <- newfiles[idx] 
  seed.num <- sapply(files, function(x) gsub(".*_seed|_cv.*", "", x))
  seed.num <- as.numeric(unname(seed.num))
  
  if(type == "rdata"){
    df <- tibble(seed = seed.num) %>%
      mutate(df = map(files,
                      ~ read_rdata(file.path(Path2Results, .), objName))) %>% 
      unnest(df)
  } else {
    df <- tibble(seed = seed.num) %>%
      mutate(df = map(files,
#                      ~ read_csv(file.path(Path2Results, .)))) %>%
                       ~ data.table::fread(file = file.path(Path2Results, .), stringsAsFactors = F, header=T ))) %>% 
      unnest(df)
  }
  
  ## some jobs might fail, so remove those
  if( 'cv' %in% colnames(df) && requireCPrevalidation) {
  num.cv <- length(unique(df$cv))
  success_seeds <- df %>% 
    group_by(seed) %>% 
    summarise(num_cv = length(unique(cv))) %>% 
    filter(num_cv == num.cv) %>% 
    pull(seed)
  df.success <- df %>% 
    filter(seed %in% success_seeds)
  
    return(df.success)
  }else{
    return(df)
  }
}

#' caret's confusionMatrix for getting precision, recall, F1, etc.
#' @param pred factor level of predictions
#' @param predprob float number for the prediction probabilities (for class 0)
#' @param ref factor level of true labels
#' 
#' @importFrom caret confusionMatrix

confm <- function(pred, predprob, ref){
  
  lvs <- levels(pred)
  
  if(length(lvs)==2){
    
    lv1 <- confusionMatrix(pred, ref, positive = lvs[1])
    lv2 <- confusionMatrix(pred, ref, positive = lvs[2])
    acc <- unname(lv1$overall[1]) ## extract accuracy
    cf <- rbind(lv1$byClass, lv2$byClass)
    rownames(cf) = paste0("Class: ", lvs)
    cf <- data.frame(cf)
    
    ## get acc and auc (auc is only for binary)
    cf["Accuracy"] <-  acc
    cf["ROCAUC"] <- pROC::auc(roc(ref, predprob))
    
  } else {
    
    cf <- confusionMatrix(pred, ref)
    acc <- unname(cf$overall[1]) ## extract accuracy
    cf <- data.frame(cf$byClass)
    cf["Accuracy"] <- acc
    
  }
  
  cf$Classes <- rownames(cf)
  
  return(cf)
}

#' Calculate root-mean-square-error, R^2, and mean-absolute-error for regression prediction
#' It also calculates the rmse if prediction is complete random
#' This is done by shuffling prediction values and recalculate rmse again
#' 
#' @param pred numerical continous values from the prediction
#' @param ref  numerical continous values from true outcome

rmse <- function(pred, ref){
  rmse <- caret::postResample(pred, ref)
  rmse <- data.frame(t(rmse))
  
  ## random shuffle score
  rmse_rnd <- caret::postResample(sample(pred), ref)
  rmse_rnd <- data.frame(t(rmse_rnd))
  colnames(rmse_rnd) <- paste0("Random_", colnames(rmse_rnd))
  
  rmse_final <- cbind(rmse, rmse_rnd)
  return(rmse_final)
}

myCindex <- function(...){
 v <- parseMajorMinor ( installed.packages()['randomForestSRC', 'Version'] )
 stopifnot(v['major'] == 2)
 if( v['minor'] >= 9 ) {
   randomForestSRC::get.cindex(...)
 }else{ # for certain earlier version of randomForestSRC
   randomForestSRC::cindex(...)
 }
}


#' Calculate hazard ratio using coxph and cindex using rfsrc::cindex
#' @param df agrregated prevalidation dataframe from getResults function
#' 
coxHR <- function(df){
  coxHR <- coxph(Surv(col_surv, col_event) ~ pred.binary, df) %>% 
    tidy %>% 
    select(estimate, p.value, conf.low, conf.high) %>%
    mutate_at(vars(estimate, conf.low, conf.high), exp) %>% 
    rename(HR = estimate)
  
  coxHR$cindex <- myCindex(df$col_surv, df$col_event, df$pred*-1)
  return(coxHR)
}

#' Calculate either prevalidation HR scores or all HR scores in all CVs
#' @param df agrregated prevalidation dataframe from getResults function
#' @param prevalid whether to use prevalidation or compute score for each CV

getHR <- function(df, prevalid = TRUE){
  if(prevalid){
    df %>% 
      group_by(seed, size) %>% 
      mutate(pred.binary = ifelse(pred < median(pred), 1, 0)) %>% 
      do(coxHR(.)) %>% 
      ungroup()
  } else {
    df %>% 
      group_by(seed, cv, size) %>% 
      mutate(pred.binary = ifelse(pred < median(pred), 1, 0)) %>% 
      do(coxHR(.)) %>% 
      ungroup()
  }
}

#' Output calculated root-mean-square-error, R^2, and mean-absolute-error 
#' on to the group_by df.preval dataframe.
#'
#' @param df prevalidation dataframe from getResults function
#' @param label_name column name of the true label
#' @param prevalid whether to use prevalidation or compute score for each CV

getRMSE <- function(df, label_name, prevalid = TRUE){
  if(prevalid){
    df %>% 
      group_by(seed, size) %>% 
      do(rmse(.$pred, .[[label_name]])) %>% 
      ungroup()
  } else {
    df %>% 
      group_by(seed, cv, size) %>% 
      do(rmse(.$pred, .[[label_name]])) %>% 
      ungroup()
  }
}

#' Calculate either prevalidation scores all scores in all CVs
#' @param df agrregated prevalidation dataframe from getResults function
#' @param label_name column name of the true label
#' @param prevalid whether to use prevalidation or compute score for each CV
#' @param metrics confusionMatrix's byClass metrics. Default: Precision, Recall, F1
#' @param calibration indicate whether to use calibrated result or not

getScores <- function(df, label_name, prevalid = TRUE, 
                      metrics = c("Precision", "Recall", "F1", "ROCAUC", "Accuracy"),
                      calibration = FALSE){
  
  ## check if not binary, remove "ROCAUC"
  if(length(unique(df[[label_name]])) > 2){
    metrics = metrics[metrics != "ROCAUC"]
  }
  
  alllevels <- levels(as.factor(df[[label_name]]))
  lowerlevel <- alllevels[1]
  
  ## use calibrated prediction 
  if(calibration){
    pred <- "pred_cali"
    predprob <- paste0("predprob_cali_", lowerlevel)
  } else {
    pred <- "pred"
    predprob <- paste0("predprob_", lowerlevel)
  }
  
  if(prevalid){
    df.scores <- df %>% 
      group_by(seed, size) %>% 
      do(confm(factor(.[[pred]], levels = alllevels), 
               .[[predprob]], factor(.[[label_name]], 
                                     levels = alllevels))) %>% 
      select(seed, size, Classes, one_of(metrics)) %>% 
      ungroup()
  } else {
    df.scores <- df %>% 
      group_by(seed, cv, size) %>% 
      do(confm(factor(.[[pred]], levels = alllevels), 
               .[[predprob]], factor(.[[label_name]], 
                                     levels = alllevels))) %>% 
      select(seed, cv, size, Classes, one_of(metrics)) %>% 
      ungroup()
  }
  
  return(df.scores)
}

#' Plotting hazard ratio barplots
#' 
#' @param df dataframe generated by getResults function
#' @param prevalid indicate whether to perform prevalidation (avg scores across CVs per seed)
#' 
#' @import ggplot2
#' 
plotHR <- function(df, prevalid = TRUE, Brier = FALSE){
  df.scores <- getHR(df, prevalid)
  metrics2show <- c("HR", "cindex")
  
  if(Brier){
    source("~/EVE/eve/reports/utils/getBrierScore.R")
    metrics2show <- c(metrics2show, "BrierScore")
    if(prevalid){
      df.brier <- get_Brier(df)
      df.scores <- df.scores %>% 
        left_join(df.brier, by=c("seed", "size"))
    } else {
      stopifnot("ibs" %in% colnames(df))
      
      df.brier <- df %>% 
        select(seed, size, cv, ibs)
      df.scores <- df.scores %>% 
        left_join(df.brier, by=c("seed", "size", "cv"))
    }
  }
  
  if(prevalid){
    plt1 <- df.scores %>% 
      select(seed, size, one_of(metrics2show)) %>% 
      gather(metrics, score, -one_of(c("seed", "size"))) %>% 
      ggplot(aes(x=as.factor(size), y=score)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("Prevalidation score across seeds") +
      xlab("# feature retain") +
      ggtitle('Prevalidation scores during RFE') +
      facet_wrap(.~metrics, dir="v", scales = "free")
  } else {
    plt1 <- df.scores %>% 
      select(seed, cv, size, one_of(metrics2show)) %>% 
      gather(metrics, score, -one_of(c("seed", "cv", "size"))) %>% 
      ggplot(aes(x=as.factor(size), y=score)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("All scores across CVs in all seeds") +
      xlab("# feature retain") +
      ggtitle('All scores during RFE') +
      facet_wrap(.~metrics, dir="v", scales = "free")
  }
  
  return(list("overall" = plt1,
              "df.scores" = df.scores))
}

plotRMSE <-  function(df, label_name, prevalid = TRUE){
  
  df.scores <- getRMSE(df, label_name, prevalid = prevalid)
  
  if(prevalid){
    ## for each class
    plt1 <- df.scores %>% 
      #select(seed, size, one_of(c("RMSE", "MAE", "Rsquared"))) %>% 
      #select_if(grepl("Random", names(.))) %>% 
      gather(metrics, score, -one_of(c("seed", "size"))) %>% 
      mutate(Rand = ifelse(grepl("Random", metrics), "Rand", "Pred")) %>% 
      mutate(metrics = gsub("Random_", "", metrics)) %>% 
      ggplot(aes(x=as.factor(size), y=score, fill = Rand)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("Prevalidation score across seeds") +
      xlab("# feature retain") +
      ggtitle('Prevalidation scores during RFE') +
      facet_wrap(.~metrics, dir="v", scales = "free")
    
  } else {
    ## for each class
    plt1 <- df.scores %>% 
      #select(seed, size, cv, one_of(c("RMSE", "MAE", "Rsquared"))) %>% 
      #select_if(grepl("Random", names(.))) %>%  
      gather(metrics, score, -one_of(c("seed", "cv", "size"))) %>% 
      mutate(Rand = ifelse(grepl("Random", metrics), "Rand", "Pred")) %>% 
      mutate(metrics = gsub("Random_", "", metrics)) %>%
      ggplot(aes(x=as.factor(size), y=score, fill = Rand)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("All scores across CVs in all seeds") +
      xlab("# feature retain") +
      ggtitle('All scores during RFE') +
      facet_wrap(.~metrics, dir="v", scales = "free")
  }
  
  return(list("overall" = plt1,
              "df.scores" = df.scores))
}


#' Plotting ML model scores barplots
#' 
#' @param df dataframe generated by getResults function
#' @param label_name column name of the true label
#' @param prevalid indicate whether to average scores across CVs per seed
#' @param metrics confusionMatrix's byClass metrics. Default: Precision, Recall, F1
#' @param calibration indicate whether to use calibrated result or not
#' 
#' @import ggplot2
#' 
plotScores <- function(df, label_name, prevalid = TRUE,
                       metrics = c("Precision", "Recall", "F1", "ROCAUC", "Accuracy"),
                       calibration = FALSE){
  
  ## check if not binary, remove "ROCAUC"
  if(length(unique(df[[label_name]])) > 2){
    metrics = metrics[metrics != "ROCAUC"]
  }
  
  df.scores <- getScores(df, label_name, prevalid = prevalid, metrics = metrics, calibration)
  
  if(prevalid){
    ## for each class
    plt1 <- df.scores %>% 
      select(seed, size, Classes, one_of(c("Precision", "Recall", "F1"))) %>% 
      gather(metrics, score, -one_of(c("seed", "size", "Classes"))) %>% 
      ggplot(aes(x=as.factor(size), y=score, color=Classes)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("Prevalidation score across seeds") +
      xlab("# feature retain") +
      ggtitle('Prevalidation scores during RFE') +
      facet_wrap(.~metrics, dir="v")
    
    ## for overall performance
    plt2 <- df.scores %>% 
      select(seed, size, Classes, one_of(c("ROCAUC", "Accuracy"))) %>% 
      filter(Classes == unique(df.scores$Classes)[1]) %>% 
      gather(metrics, score, -one_of(c("seed", "size", "Classes"))) %>% 
      ggplot(aes(x=as.factor(size), y=score)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("Prevalidation score across seeds") +
      xlab("# feature retain") +
      ggtitle('Prevalidation scores during RFE') +
      facet_wrap(.~metrics, dir="v")
    
  } else {
    ## for each class
    plt1 <- df.scores %>% 
      select(seed, cv, size, Classes, one_of(c("Precision", "Recall", "F1"))) %>% 
      gather(metrics, score, -one_of(c("seed", "cv", "size", "Classes"))) %>% 
      ggplot(aes(x=as.factor(size), y=score, color=Classes)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("All scores across CVs in all seeds") +
      xlab("# feature retain") +
      ggtitle('All scores during RFE') +
      facet_wrap(.~metrics, dir="v")
    
    ## for overall performance
    plt2 <- df.scores %>% 
      select(seed, cv, size, Classes, one_of(c("ROCAUC", "Accuracy"))) %>% 
      gather(metrics, score, -one_of(c("seed", "cv", "size", "Classes"))) %>% 
      ggplot(aes(x=as.factor(size), y=score)) + 
      geom_boxplot() +
      theme_bw() +
      ylab("All scores across CVs in all seeds") +
      xlab("# feature retain") +
      ggtitle('All scores during RFE') +
      facet_wrap(.~metrics, dir="v")
    
  }
  
  return(list("byClass" = plt1,
              "overall" = plt2,
              "df.scores" = df.scores))
}

#' Summary table of min/max median of each performance metrics
#' @param df data.frame of df.scores from getScores function
#'  
eval.minmax <- function(df){
  cols2rm <- c("seed", "size", "p.value", "conf.low", "conf.high", "Classes",
               "Random_MAE", "Random_RMSE", "Random_Rsquared")
  
  score.max <- df %>% 
    gather("metrics", "score", -one_of(cols2rm)) %>% 
    group_by(seed, size, metrics) %>% 
    summarize(avg.score = mean(score, na.rm = T)) %>% ## average across CV (if not prevalidation)
    group_by(size, metrics) %>% 
    summarize(med = median(avg.score, na.rm = T)) %>% 
    group_by(metrics) %>% 
    arrange(desc(med)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    rename(median.max = med,
           size.max = size) 
  
  score.min <- df %>% 
    gather("metrics", "score", -one_of(cols2rm)) %>% 
    group_by(seed, size, metrics) %>% 
    summarize(avg.score = mean(score, na.rm = T)) %>% ## average across CV (if not prevalidation)
    group_by(size, metrics) %>% 
    summarize(med = median(avg.score, na.rm = T)) %>% 
    group_by(metrics) %>% 
    arrange(med) %>% 
    slice(1) %>% 
    ungroup() %>% 
    rename(median.min = med,
           size.min = size)
  scores <- score.max %>% 
    left_join(score.min, by="metrics") %>% 
    select(metrics, size.max, median.max, size.min, median.min)
  
  return(scores)
}

#' Plot calibration curve based on given feature size
#' This plot is only for binary classification, although multiclass can be classified in sklearn.
#' 
#' @param df agrregated prevalidation dataframe from getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' 
plotCalibration <- function(df, ft_num = NULL, cuts = 10){
  
  if(length(unique(df[[label_name]])) > 2){
    stop("Currently the plot only works for binary class problem.")
  }
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  }
  
  lowerlevel <- levels(as.factor(df$pred))[1]
  preprob_cali_name <- paste0("predprob_cali_", lowerlevel)
  preprob_name <- paste0("predprob_", lowerlevel)
  
  df.in <- df %>% filter(size == ft_num)
  cal_obj <- calibration(as.factor(eval(as.name(label_name))) ~ eval(as.name(preprob_cali_name)) + eval(as.name(preprob_name)),
                         data = df.in, cuts = cuts)
  plot(cal_obj, type = "l", auto.key = list(columns = 3,
                                            lines = TRUE,
                                            points = FALSE))
}

#' Plot VIMP for non-XGBoost models
#' This is based on 'vimp' column of the input data
#' 
#' @param df dataframe generated by getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' @param top_n integer, choose number of top VIMP to show. Default is top 20.
#' @param top_n_by size or freq. Select top_n features according to average coefficient size, or according to frequency of occurrence across CVs(before averging within a seed). Frequency may be preferred if X was not scaled.
#' @param ft_name character vector, to limit VIMP plot to certain features. When ft_name is provided, top_n will be ignored
#' @param bins integer, number of bins for histogram of feature importance frequency
#' @param ignore_neg bool, choose whether to use the original value of vimp when assessing variable importance 
#'                   It operates on average (across seed) of vimp sum(across CV within seed).
#'                   In lasso, negative coefficient is important as well (ignore_neg= F, default). 
#'                   But in rfsrc, negative vimp might be least important (ignore_neg=T) and hence ignored. 
#' @param non_zero_value a small positive real number (eg 1e-6). If this value is specified, records with abs(vimp) less than the value will be removed before any calculation. This is mainly to support multi-variate lasso where certain features with zero coefficients are included.
#' 
plotVIMP2 <- function(df, ft_num=NULL, top_n=20, top_n_by="freq", ft_name=NULL,
                      bins=50, ignore_neg=FALSE, non_zero_value=NA){
  stopifnot(top_n_by %in% c("size", "freq"))
  if(nrow(df) == 0) {
    cat('no feature retained')
    return(NULL)
  }

  if(is.null(ft_num)){
    ft_num <- max(df$size)
  }

  # only focus on selected size
  # also needed for `rank` calculations below
  # otherwise, get error like:
  # Error: Column `rank` must be length 15806 (the number of rows) or one, not 15801
  # > traceback()
  # 2: df %>% count(feature) %>% arrange(desc(n)) %>% mutate(rank = 1:nrow(df.vimp.scores)) at PerformanceUtils.R#699
  # 1: plotVIMP2(df.vimp, bin = 100, ignore_neg = T, ft_num = ft_num)
  df <- df %>%
    filter(size == ft_num)
  
  if(! (missing(non_zero_value) || is.null(non_zero_value) ||
        is.na(non_zero_value) || non_zero_value <= 0)){
    toRemove <- abs( df[, 'vimp'] ) < non_zero_value
    print(paste('removing', sum(toRemove),'records as their vimp is less than',
                non_zero_value))
    df <- df[!toRemove, ,drop=F]
  }

  df.vimp.scores <- df %>% 
    filter(size == ft_num) %>% 
    group_by(seed, feature) %>% 
    summarise(sum_vimp = sum(vimp, na.rm = T)) %>% 
    group_by(feature) %>% 
      summarise(
          ## derive mean and CI for plots
          vimp.avg = mean(sum_vimp, na.rm = T),
          lower = quantile(sum_vimp, probs=0.025),
          upper = quantile(sum_vimp, probs=0.975)
          ) %>% 
    ## sign is for lasso's coefficient which may contain negative values
    mutate(Sign = ifelse(vimp.avg>=0, "Pos", "Neg"), 
           vimp.avg.abs = abs(vimp.avg))
  
  if(ignore_neg){
    df.vimp.scores <- df.vimp.scores %>% 
      mutate(vimp.avg.abs = vimp.avg)
  }
  
  ## used if top_n_by="freq"
  ## start with ORIGINAL df
  df.feature.freq <- df %>%
      ## add column `n` with counts
      count(feature) %>%
      arrange(desc(n)) %>%
      ## also derive rank if preferred for plots
      ## nrow(df.vimp.scores) == number of features
      mutate(rank=1:nrow(df.vimp.scores))

  assertthat::assert_that(nrow(df.feature.freq) == nrow(df.vimp.scores))

  ## User can provide which features they want in `ft_name`. Otherwise, take
  ## `top_n` features according to `top_n_by`
  df.top.f <- df.vimp.scores %>% 
      ## add rank info for top_n filtering and plotting
      left_join(
          df.feature.freq %>% select(feature, n, rank),
          by="feature")
  
  if(!is.null(ft_name)) {
    df.top.f <- subset(df.top.f, feature %in% ft_name )
  } else if (top_n_by == "size") {
    df.top.f <- df.top.f %>% 
      # also arrange by most frequent if tied
      arrange(desc(vimp.avg.abs), desc(n)) %>%
      slice(1:top_n) 
  } else if (top_n_by == "freq") {
      df.top.f <- df.top.f %>%
          ## for ridge, all freq's are the same, so also arrange by size
          arrange(desc(n), desc(vimp.avg.abs)) %>%
          slice(1:top_n)
  }

  ## PLOT 1 of 2
  plt.dist.f <- ggplot(df.vimp.scores, aes(x=vimp.avg)) + 
    geom_histogram(color="black", fill="white", bins=bins) +
    theme_bw() +
    ylab("Frequency") +
    xlab("Average (across seeds) of vimp (sum of all CVs per seed)") +
    ggtitle(paste("From", ft_num, "feature step based on vimp"))
  
  ## PLOT 2 of 2
  ## extract options for plot annotation
  RANGE <- diff(range(df.top.f$vimp.avg, finite = TRUE))
  YMIN <- min(df.top.f$vimp.avg) - RANGE/2
  ## for plotting, order features according to data frame order (as arrange()'d
  ## above)
  df.top.f.ggplot <- df.top.f %>%
      mutate(feature = factor(feature, levels=feature)) %>%
      ## reverse factor order to match with coord_flip() below
      ## https://stackoverflow.com/q/34227967/3217870
      mutate(feature = forcats::fct_rev(feature))

  if (top_n_by == "size") {
      plt.features.base <- ggplot(df.top.f.ggplot, aes(x=feature, y=vimp.avg)) +
          ggtitle(paste("Feature importance by",
                        ifelse(ignore_neg,'','absolute')," VIMP magnitude"))
  }
  else {
      plt.features.base <- ggplot(df.top.f.ggplot, aes(x=feature, y=vimp.avg)) +
          ggtitle("Feature importance by usage frequency")
  }

  plt.features <- plt.features.base +
      geom_col() +
      geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, alpha=1) +
      geom_text(aes(y = YMIN, label=n)) +
      labs(x="Feature", y="Importance (95% CI)") +
      ## this flips the order, so using forcats::fct_rev() above
      coord_flip() +
      theme_bw()
      
  return(list(df = df.vimp.scores,
              plt.dist.f = plt.dist.f,
              plt.features = plt.features,
              top_df = df.top.f
              ))
  }
  


#' Plot Coefficient for lasso. Modified from plotVIMP2
#' This is based on 'coef' column of the input data
#' 
#' @param df dataframe generated by getResults function
#' @param top_n integer, choose number of top coefficients to show. Default is top 20.
#' @param top_n_by coef or freq. Select top_n features according to average coefficient, or according to frequency of occurrence across CVs (before aggregating within a seed). Frequency may be preferred if X was not scaled.
#' @param ft_name character vector, to limit the plot to certain features. When ft_name is provided, top_n and top_n_by will be ignored
#' @param ci1 1 sided confidence level for error bars
#' 
plotCoef <- function(df, top_n=20, top_n_by="freq", ft_name=NULL, ci1=0.025){
  
# modified from plotVIMP2. Removed the previous argument of ft_num, ignore_neg and non_zero_value. These just makes the function hard to read

  #' ignore_neg bool, choose whether to use the original value of vimp when assessing variable importance. It operates on average (across seed) of vimp sum(across CV within seed).
  #' In lasso, negative coefficient is important as well (ignore_neg= F, default). 
  #'But in rfsrc, negative vimp might be least important (ignore_neg=T) and hence ignored.     


  stopifnot( all(c('coef') %in% colnames(df) ) )
  
  if(nrow(df) == 0) {
    cat('no feature retained')
    return(NULL)
  }
  
  if('size' %in% colnames(df)) stopifnot(length(unique(df$size)) ==1)

  sum_coef <-  with(df, tapply(coef, list(seed, feature), sum)) 
  coef.avg <- colMeans(sum_coef,na.rm=T)
  df.coef.score <- data.frame(coef.avg = coef.avg, Sign = ifelse(coef.avg>=0, "Pos", "Neg"), 
                              coef.avg.abs = abs(coef.avg),
                              feature = names(coef.avg), 
                              lower = apply(sum_coef, 2, quantile, probs=ci1, na.rm=T),
                              upper = apply(sum_coef, 2, quantile, probs= 1- ci1, na.rm=T), stringsAsFactors = F)
  
   
  df.feature.freq <- data.frame( table(df$feature) )
  colnames(df.feature.freq) <- c('feature','n')
  df.feature.freq$feature <- as.character(df.feature.freq$feature)
#  df.feature.freq <- df.feature.freq[order(df.feature.freq$n * -1), ]
#  df.feature.freq$rank_n <- 1:nrow(df.feature.freq)
  
  
  if( 'vimp' %in% colnames(df)){
    sel_vimp_na <- is.na(df$vimp)
    if(any(sel_vimp_na)){
      min_vimp <- min(df$vimp[!sel_vimp_na] )
      print(paste('there are', sum(sel_vimp_na), 'NA values in vimp before summation within seeds; they are imputed with the smallest vimp value of',  ifelse(min_vimp >0, 0, round.signif(min_vimp,2))))
      df[sel_vimp_na,'vimp'] <- min(0, min_vimp)
      rm(sel_vimp_na, min_vimp)
    }
    
    sum_vimp <-   with(df, tapply(vimp, list(seed, feature), sum))
  
    sel_vimp_na <- is.na(sum_vimp)
    if( any(as.vector(sel_vimp_na)) ){
      min_vimp <- min(as.vector(sum_vimp[!sel_vimp_na]) )
      print(paste('there are', sum(as.vector(sel_vimp_na)), 'NA values in vimp after summation within seeds; they are imputed with the smallest vimp value of',  ifelse(min_vimp >0, 0, round.signif(min_vimp,2))))
      sum_vimp[sel_vimp_na] <- min_vimp
    }
  
  vimp.avg <- colMeans(sum_vimp)
 stopifnot(!any(is.na(vimp.avg)))
 
  df.feature.freq$vimp <- vimp.avg[as.character(df.feature.freq$feature)]
  }
#  df.feature.freq$'vimp.avg' <- df.feature.freq$vimp
  
#  df.feature.freq <- df.feature.freq[order(df.feature.freq$vimp * -1), ]
#  df.feature.freq$rank_vimp <- 1:nrow(df.feature.freq)
  
  
  assertthat::assert_that(nrow(df.feature.freq) == nrow(df.coef.score))
  assertthat::assert_that(all(sort(df.feature.freq$feature) == sort(df.coef.score$feature)))
  
  
  df <- df.top.f <- merge(df.coef.score,  df.feature.freq, by='feature')
  
  stopifnot( top_n_by %in% c("coef", "freq", "gain"))
  
  plt.dist <- NULL     ## PLOT 1 of 2
  
  if(!is.null(ft_name)) {
    df.top.f <- subset(df.top.f, feature %in% ft_name )
  } else {
    if (top_n_by == "coef") {
      
      plt.dist  <- ggplot(df.coef.score, aes(x=coef.avg)) + 
        #    geom_histogram(color="black", fill="white", bins=bins) +
        geom_histogram(color="black", fill="white") +
        theme_bw() +
        ylab("Frequency") +
        xlab("Average (across seeds) of coefficient (summed across all CVs per seed)") +
        ggtitle(paste("Distribution across all", nrow(df.coef.score),'features'))
      
      df.top.f <- df.top.f %>% 
        # also arrange by most frequent if tied
        arrange(desc(coef.avg.abs), desc(n) ) %>%
        slice(1:top_n) 
    } else {
      if (top_n_by == "freq") {
        
        plt.dist <- ggplot(df.feature.freq, aes(x=n)) + 
          #    geom_histogram(color="black", fill="white", bins=bins) +
          geom_histogram(color="black", fill="white") +
          theme_bw() +
          ylab("Frequency") +
          xlab("Frequency of use across CVs and seeds") +
          ggtitle(paste("Distribution across all", nrow(df.coef.score),'features'))
        
        df.top.f <- df.top.f %>%
          ## for ridge, all freq's are the same, so also arrange by size
          arrange(desc(n), desc(coef.avg.abs)) %>%
          slice(1:top_n)
      } else {
        if (top_n_by == "gain") {
          stopifnot( 'vimp' %in% colnames(df))
          plt.dist <- ggplot(df.feature.freq, aes(x=vimp)) + 
            #    geom_histogram(color="black", fill="white", bins=bins) +
            geom_histogram(color="black", fill="white") +
            theme_bw() +
            ylab("Frequency") +
            xlab("Average (across seeds) of worsen statistics\nfrom NextDoor analysis (summed across all CVs per seed)") +
            ggtitle(paste("Distribution across all", nrow(df.coef.score),'features'))
          
      df.top.f <- df.top.f %>%
        arrange(desc(vimp), desc(n)) %>%
        slice(1:top_n)
        }
      }
      }
  }
  
  ## PLOT 2 of 2
  ## extract options for plot annotation
  RANGE <- diff(range(df.top.f$coef.avg, finite = TRUE))
  YMIN <- min(df.top.f$coef.avg) - RANGE/2
  ## for plotting, order features according to data frame order (as arrange()'d
  ## above)
  df.top.f.ggplot <- df.top.f %>%
    mutate(feature = factor(feature, levels=feature)) %>%
    ## reverse factor order to match with coord_flip() below
    ## https://stackoverflow.com/q/34227967/3217870
    mutate(feature = forcats::fct_rev(feature))
  
  plt.features.base <- ggplot(df.top.f.ggplot, aes(x=feature, y=coef.avg))
  
  if (top_n_by == "coef") {
    plt.features.base <- plt.features.base +
      ggtitle(paste("Top feature, by absolute Coeficient magnitude"))
  }  else {
    if (top_n_by == "freq") {
    plt.features.base <- plt.features.base +
      ggtitle("Top feature, by usage frequency")
    }else{
      if (top_n_by == "gain") {
      plt.features.base <- plt.features.base +
        ggtitle("Top feature, by the worsen statistic from NextDoor analysis")
      }else{
        plt.features.base <- plt.features.base +
          ggtitle("Top features")
      }
    }}
  
  plt.features <- plt.features.base +
    geom_col() +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4, alpha=1) +
    geom_text(aes(y = YMIN, label=n)) +
    labs(x="Feature", y=paste("Coefficent (", round( (1 - ci1 * 2) *100, 0) ,"% CI)")) +
    ## this flips the order, so using forcats::fct_rev() above
    coord_flip() +
    theme_bw()
  
  return(list(df =df , 
              plt.dist  = plt.dist, plt.features = plt.features,
              top_df = df.top.f
  ))
}



#' This is to visualize Lasso's multi-class prediciton VIMP
#' Because it outputs vimp for each class, we have to loop through each class
#' 
plotVIMP2.multiclass <- function(df, ft_num=NULL, top_n=20, top_n_by="freq",
                                 bins=30, ...){
  idx <- grepl("vimp_", colnames(df))
  vimps <- colnames(df)[idx]
  classes <- gsub("vimp_", "", vimps)
  
  for(i in 1:length(vimps)){
    df.vimp.tmp <- df %>% 
      select(-one_of(vimps[-i])) 
    idx.col <- grepl("vimp_", colnames(df.vimp.tmp))
    colnames(df.vimp.tmp)[idx.col] <- "vimp"
    
    df.vimp.plt <- plotVIMP2(df.vimp.tmp, bins = bins, top_n = top_n,
                             top_n_by = top_n_by, ...)
    ## plot
    plt1 <- df.vimp.plt$plt.dist.f +
      ggtitle(paste("Distribution of Feature Coefficient for Class", classes[i]))
    print(plt1)
    plt2 <- df.vimp.plt$plt.features +
      ggtitle(paste("Top Features for Class", classes[i]))
    print(plt2)
  }
}

#' Plot top VIMP (vimp scores are averaged across seeds)
#' This is based on 'frequency' of the feature existing in all the trees,
#' and normalized by 'total number of nodes(features)'
#' 
#' @param df dataframe generated by getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' @param top_n integer, choose number of top VIMP to show. Default is top 20.
#' @param bins integer, number of bins for histogram of vimp frequencies
#' 
plotVIMP <- function(df, ft_num=NULL, top_n=20, bins=50){
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  }
  
  df.vimp.scores <- df %>% 
    filter(size == ft_num) %>% 
    group_by(seed, feature) %>% 
    summarise(sum_weight = sum(weight, na.rm = T),
              sum_gain   = sum(gain  , na.rm = T)) %>% 
    group_by(feature) %>% 
    summarise(freq.avg = mean(sum_weight, na.rm = T),
              gain.avg = mean(sum_gain  , na.rm = T)) 
  
  df.top.f <- df.vimp.scores %>% 
    arrange(desc(freq.avg)) %>% 
    slice(1:top_n) %>% 
    arrange(freq.avg) ## this is to plot the vimp from top-down
  
  df.top.g <- df.vimp.scores %>% 
    arrange(desc(gain.avg)) %>% 
    slice(1:top_n) %>% 
    arrange(gain.avg)
  
  plt1 <- ggplot(df.vimp.scores, aes(x=freq.avg)) + 
    geom_histogram(color="black", fill="white", bins=bins) +
    theme_bw() +
    ylab("Frequency") +
    xlab("average (across seeds) of Frequency (sum of all CVs per seed)") +
    ggtitle(paste("with", ft_num, "features based on Frequency"))
  
  plt2 <- ggplot(df.vimp.scores, aes(x=gain.avg)) + 
    geom_histogram(color="black", fill="white", bins=bins) +
    theme_bw() +
    ylab("Frequency") +
    xlab("average (across seeds) of Gain (sum of all CVs per seed)") +
    ggtitle(paste("with", ft_num, "features based on Gain"))
  
  plt3 <- ggplot(df.top.f, aes(x=feature, y=freq.avg)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits=df.top.f$feature) +
    coord_flip() +
    theme_bw() +
    ylab("Feature Importance (frequency)") +
    ggtitle(paste("Top", top_n, "features at", ft_num, "feature set based on Frequency"))
  
  plt4 <- ggplot(df.top.g, aes(x=feature, y=gain.avg)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits=df.top.g$feature) +
    coord_flip() +
    theme_bw() +
    ylab("Feature Importance (gain)") +
    ggtitle(paste("Top", top_n, "features at", ft_num, "feature set based on Gain"))
  
  return(list(df = df.vimp.scores,
              plt.dist.f = plt1,
              plt.dist.g = plt2,
              plt.fts.f  = plt3,
              plt.fts.g  = plt4))
}

#' Plot top VIMP (Gain vs Frequency)
#' 
#' These scores are averaged across seeds.
#' 'Frequency': the feature existing in all the trees (weight is the variable used in xgboost)
#' 'Gain': the amount of increase in score when the feature is added to split
#' 
#' @param df dataframe generated by getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' @param top_n integer, choose number of top N Gain to show. Default is top 20.
#' 
#' @import ggrepl
#' 
plotVIMP_scatter <- function(df, ft_num=NULL, top_n=20){
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  }
  
  df.vimp.topN <- df %>% 
    filter(size == ft_num) %>% 
    group_by(seed, feature) %>% 
    summarise(sum_weight = sum(weight, na.rm = T),
              sum_gain   = sum(gain  , na.rm = T)) %>% 
    group_by(feature) %>% 
    summarise(freq.avg = mean(sum_weight, na.rm = T),
              gain.avg = mean(sum_gain  , na.rm = T)) %>% 
    arrange(desc(gain.avg), desc(freq.avg)) %>%
    #mutate_at(vars(gain.avg), funs(as.numeric(scale(.)))) %>%  ## rescale the data (0 mean/var)
    slice(1:top_n)
  
  ggplot(df.vimp.topN, aes(freq.avg, gain.avg, label = feature)) +
    geom_point(alpha = 0.7, shape = 1) +
    ggrepel::geom_text_repel(
      nudge_x = 0,
      nudge_y = 0,
      segment.alpha = 0.5) +
    theme_bw() +
    xlab("Frequency (averaged across seeds)") +
    ylab("Gain (averaged across seeds)") +
    ggtitle(paste("Top", top_n, "features at", ft_num, "feature set"))
  
}

#' Plot grid search results (f1 scores are averaged across seeds)
#' 
#' @param df dataframe generated by getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' @param max_seed_num integer, total number of plots to visualize, each plot represent one seed
#' 
plotGridS <- function(df, ft_num=NULL, max_seed_num=5){
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  } else{
    ft_num <- ft_num
  }

  seeds <- unique(df.grid$seed)
  
  if(length(seeds)>max_seed_num){
    seeds <- seeds[1:max_seed_num]
  }
  
  plt <- df %>% 
    filter(size == ft_num,
           seed %in% seeds) %>%
    mutate_at(vars(-one_of(c("seed", "size", "cv", "score"))),
              funs(scales::rescale)) %>% 
    gather(hyperP, values, -one_of(c("seed", "size", "cv", "score"))) %>%
    ggplot(aes(x=hyperP, y=values)) + 
    geom_line(aes(color = score, group = cv),
              size=1, alpha=0.8) +
    theme_bw() +
    facet_wrap(~seed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_viridis_c() +
    ylab("Normalized Hyperparamter Values ") +
    xlab("Hyperparameters") +
    ggtitle('Hyperparamter Tuning in each CV',
            subtitle = paste0("each box represents different seed (maximum ", max_seed_num, " seeds are shown)"))
  
  return(plt)
}

#' Generate Confustion Matrix for classification task
#' @param df dataframe generated by getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' 
confusionMat <- function(df, ft_num=NULL){
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  }
  num_seed <- length(unique(df$seed))
  
  df.tmp <- df %>% 
    filter(size == ft_num)
  
  cm <- confusionMatrix(as.factor(df.tmp$pred), 
                        as.factor(df.tmp[[runSpec$label_name]]))
  cat("confusion matrix at feature size =", ft_num, "\nsum across", num_seed, "seeds\n")
  cm$table
}


correlation <- function(df, seed_num=NULL){
  
  if(is.null(seed_num)){
    seed_num <- min(df$seed)
  }
  
  df.in <- df %>% filter(seed == seed_num)
  df.in <- subset(df.in, size == max(size))
  lb_name <- as.name(runSpec$label_name)
  plt <- ggplot(data = df.in, 
                aes_string(x="pred", y = lb_name)) + 
    geom_point(alpha=0.3) +
    theme_bw() +
    ggtitle(paste0("Correlation at seed = ", seed_num,' using ', max(df.in$size),' feature set input'))
  
  lb_name <- as.name(runSpec$label_name)
  
  cor.summary <- df %>% 
    group_by(seed) %>% 
    summarise(corr = cor(!!lb_name, pred)) %>% 
    summarise(cor.avg = mean(corr),
              cor.sdt = sd(corr))
  
  print(knitr::kable(cor.summary, caption="Averaged pearson correlation across seeds"))
  print(plt)
}

#' Plot feature importance heatmap
#'  @param df dataframe of average vimp. usually from the output of plotVIMP/plotVIMP2
#'  @param df.orig dataframe of original input data
#' 
plotHeatmap <- function(df, df.orig, top_n=20, classification=FALSE, scaleBy='row'){
  
  
  if(grepl("lasso", tolower(runSpec$engine))){
#    vimp_name <- as.name("vimp.avg")
#    vimp_sort <- as.name("vimp.avg.abs")
    vimp_name <- as.name("vimp")
    vimp_sort <- as.name("vimp")
    
      } else if(grepl("xgboost", tolower(runSpec$engine))){
    vimp_name <- as.name("gain.avg")
    vimp_sort <- as.name("gain.avg")
  } else if(grepl("rfesrc", tolower(runSpec$engine))){
    vimp_name <- as.name("vimp.avg")
    vimp_sort <- as.name("vimp.avg")
  } else {
    stop("unrecognized engine")
  }
  
  fts <- df %>% 
    arrange(desc(!!vimp_sort)) %>% 
    slice(1:top_n) %>% 
    select(feature, !!vimp_name)
  
  if(is.na(runSpec$event_col) && !is.na(runSpec$label_name)){
    lb <- runSpec$label_name
  } else if(!is.na(runSpec$event_col) && is.na(runSpec$label_name)){
    lb <- runSpec$surv_col
  } else {
    stop("unrecognized label (Y)")
  }
  
  df.orig <- df.orig %>% arrange(!!as.name(lb))
  
  df.plt <- df.orig[, fts$feature]
  df.lb <- df.orig[, lb]
  
  df.in <- data.frame(t(df.plt))
  
  annot.row <- data.frame("vimp" = fts[[vimp_name]])
  rownames(annot.row) <- rownames(df.in)
  if(classification){
    annot.col <- data.frame(df.orig[, lb]) %>%  mutate_all(as.factor)
  } else {
    annot.col <- data.frame(df.orig[, lb])
  }
   
  #annot.col <- as.factor(annot.col)
  rownames(annot.col) <- colnames(df.in)
  title <- paste0("Heatmap of top ", top_n, " important features")
  plt <- pheatmap::pheatmap(df.in,
                            scale = scaleBy,
                            annotation_row = annot.row,
                            annotation_col = annot.col,
                            cluster_cols = F,
                            show_colnames = F,
                            annotation_names_row = F,
                            main = title)
  print(plt)
}




