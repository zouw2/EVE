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
#' 
getResults <- function(project_home, project_name, type, objName = NULL){
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
    df <- data_frame(seed = seed.num) %>% 
      mutate(df = map(files,
                      ~ read_rdata(file.path(Path2Results, .), objName))) %>% 
      unnest() 
  } else {
    df <- data_frame(seed = seed.num) %>% 
      mutate(df = map(files,
                      ~ read_csv(file.path(Path2Results, .)))) %>% 
      unnest() 
  }
  
  ## some jobs might fail, so remove those
  num.cv <- length(unique(df$cv))
  success_seeds <- df %>% 
    group_by(seed) %>% 
    summarise(num_cv = length(unique(cv))) %>% 
    filter(num_cv == num.cv) %>% 
    pull(seed)
  df.success <- df %>% 
    filter(seed %in% success_seeds)
  
  return(df.success)
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


#' Calculate hazard ratio using coxph and cindex using rfsrc::cindex
#' @param df agrregated prevalidation dataframe from getResults function
#' 
coxHR <- function(df){
  coxHR <- coxph(Surv(col_surv, col_event) ~ pred.binary, df) %>% 
    tidy %>% 
    select(estimate, p.value, conf.low, conf.high) %>%
    mutate_at(vars(estimate, conf.low, conf.high), exp) %>% 
    rename(HR = estimate)
  
  coxHR$cindex <- randomForestSRC::cindex(df$col_surv, 
                                          df$col_event, 
                                          df$pred*-1)
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
      ggplot(aes(x=as.factor(size), y=score, fill=Classes)) + 
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
      ggplot(aes(x=as.factor(size), y=score, fill=Classes)) + 
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
#' 
plotVIMP2 <- function(df, ft_num=NULL, top_n=20, bin=30){
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  }
  
  df.vimp.scores <- df %>% 
    filter(size == ft_num) %>% 
    group_by(seed, feature) %>% 
    summarise(sum_vimp = sum(vimp, na.rm = T)) %>% 
    group_by(feature) %>% 
    summarise(vimp.avg = mean(sum_vimp, na.rm = T)) %>% 
    mutate(sign = ifelse(vimp.avg>=0, "pos", "neg"), ## this is for lasso's coefficient which may contain negative values
           vimp.avg.abs = abs(vimp.avg))
  
  df.top.f <- df.vimp.scores %>% 
    arrange(desc(vimp.avg.abs)) %>% 
    slice(1:top_n) %>% 
    arrange(vimp.avg.abs) ## this is to plot the vimp from top-down
  
  plt1 <- ggplot(df.vimp.scores, aes(x=vimp.avg)) + 
    geom_histogram(color="black", fill="white", bins=50) +
    theme_bw() +
    ylab("Frequency") +
    xlab("average (across seeds) of vimp (sum of all CVs per seed)") +
    ggtitle(paste("with", ft_num, "features based on vimp"))
  
  plt2 <- ggplot(df.top.f, aes(x=feature, y=vimp.avg.abs, fill=sign)) +
    geom_bar(stat="identity", alpha=0.7) +
    scale_x_discrete(limits=df.top.f$feature) +
    scale_fill_manual(values=c("#999999", "#000000")) +
    coord_flip() +
    theme_bw() +
    ylab("Feature Importance") +
    ggtitle(paste("Top", top_n, "features at", ft_num, "feature set"))
  
  return(list(df = df.vimp.scores,
              plt.dist.f = plt1,
              plt.fts.f  = plt2))
}

#' Plot top VIMP (vimp scores are averaged across seeds)
#' This is based on 'frequency' of the feature existing in all the trees,
#' and normalized by 'total number of nodes(features)'
#' 
#' @param df dataframe generated by getResults function
#' @param ft_num integer, choose which RFE round to plot. Default is using all features
#' @param top_n integer, choose number of top VIMP to show. Default is top 20.
#' 
plotVIMP <- function(df, ft_num=NULL, top_n=20, bin=30){
  
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
    geom_histogram(color="black", fill="white", bins=50) +
    theme_bw() +
    ylab("Frequency") +
    xlab("average (across seeds) of Frequency (sum of all CVs per seed)") +
    ggtitle(paste("with", ft_num, "features based on Frequency"))
  
  plt2 <- ggplot(df.vimp.scores, aes(x=gain.avg)) + 
    geom_histogram(color="black", fill="white", bins=50) +
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
    geom_text_repel(
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
#' 
plotGridS <- function(df, ft_num=NULL){
  
  if(is.null(ft_num)){
    ft_num <- max(df$size)
  } else{
    ft_num <- ft_num
  }
  
  # plt <- df %>% 
  #   filter(size == ft_num) %>% ## select only specific RFE set
  #   gather(hyperP, values, -one_of(c("seed", "size", "cv", "f1"))) %>%
  #   group_by(seed, size, hyperP, values) %>% 
  #   summarise(avg.score = mean(f1)) %>% 
  #   ggplot(aes(x=as.factor(values), y=avg.score)) + 
  #   geom_boxplot() +
  #   theme_bw() +
  #   ylab("Averaged F1 scores across seeds") +
  #   xlab("values") +
  #   ggtitle('Averaged F1 scores during RFE') +
  #   facet_wrap(.~ hyperP)
  
  plt <- df %>% 
    filter(size == ft_num) %>%
    gather(hyperP, values, -one_of(c("seed", "size", "cv", "score"))) %>%
    ggplot(aes(x=score, y=values, color=hyperP)) + 
    geom_point(size=3, alpha=0.8) +
    theme_bw() +
    ylab("Hyperparamters values ") +
    xlab("Averaged scores across CVs and seeds") +
    labs(color='Hyperparameter') +
    ggtitle('Hyperparamters vs Scores')
  
  return(plt)
}