#' This is modified from Wei's Brier score calculation
#' 1. Takes rfsrc's predicted probabilities / sample 
#' 2. Fill in missing time points. Since different CVs might give different probilitiy time points,
#'     we need to interpolate those missing time points. 
#' 3. Calculate Brier score
#' 

#' Main functions
#' require 'runSpec' in the environment

get_Brier <- function(df.preval){
  df.survprob <- df.preval %>% 
    select(seed, cv, size, surv_prob, surv_prob.times) %>% 
    group_by(seed, size, cv) %>% 
    nest() %>%
    mutate(data = map(data, ~ surv_split(.)))
  
  df.surv.all <- data.frame()
  seeds <- unique(df.survprob$seed)
  for(s in seeds){
    df.sub <- df.survprob %>% 
      filter(seed == s) 
    
    commonCols <- get_times(df.sub)
    
    df.sub <- df.sub %>% 
      mutate(data = map(data, 
                        ~ fillInMissingProbWithNeighbor(., commonCols, method = 'LinearInterpolation'))) 
    df.surv.all <- rbind(df.surv.all, df.sub)
  }
  df.surv.all <- df.surv.all %>% unnest()
  if(is.na(runSpec$sample_ID)) runSpec$sample_ID <- "RowIndex"
  df.surv.all[[runSpec$sample_ID]] <- df.preval[[runSpec$sample_ID]]
  
  seed0 <- unique(df.preval$seed)[1]
  size0 <- unique(df.preval$size)[1]
  df.y <- df.preval %>% 
    filter(seed == seed0, size == size0) %>% 
    select(col_surv, col_event, !!as.name(runSpec$sample_ID))
  aggB <- aggBrier(data.frame(df.surv.all), 
                   data.frame(df.y),
                   outcome.var = c("col_surv", "col_event"),
                   idVar = runSpec$sample_ID)
  aggB <- reshape2::melt(aggB)
  colnames(aggB) <- c("size", "seed", "BrierScore")
  return(aggB)
}


surv_split <- function(df){
  df.out <- separate(df, surv_prob, 
                     into = df$surv_prob.times[[1]], sep=",") %>% 
    select(-surv_prob.times) %>% 
    mutate_all(as.numeric)
  return(df.out)
}

extract_num <- function(cols){
  cols.out <- gsub("^X", "", cols)
  return(as.numeric(cols.out))
}

get_times <- function(df){
  col.times <- df %>% 
    mutate(times = map(data, ~ extract_num(colnames(.)))) %>% 
    unnest(times) %>% 
    arrange(times) %>% 
    pull(times) %>% 
    unique() %>% 
    as.numeric()
  
  return(col.times)
}

fillInMissingProbWithNeighbor <- function(x, cols, pattern='^X[\\d\\.]+$', method = c('previousEstimate', 'LinearInterpolation' )){
  ## ToDO: temporarily hack, need to improve the code
  if(!grepl("^X", colnames(x)[1])){
    colnames(x) <- paste0("X", colnames(x))
    cols <- paste0("X", cols)
  }
  
  extra <- setdiff(colnames(x), cols)
  if(length(extra) > 0) stop(paste('unit run matrix contains unexpected column', paste(extra, collapse = ',')))
  
  if(all(cols %in% colnames(x))) return(x[, cols, drop=F])
  
  stopifnot(! 'one' %in% colnames(x))
  
  colN <- grep(pattern, cols, perl=T, value=T)
  
  tofill <- setdiff(colN, colnames(x))
  #  toaddBack <- setdiff(colnames(x), cols)
  if(method == 'LinearInterpolation'){
    fillin <- sapply(tofill, function(v, timeAvailable = as.numeric( gsub('^X','', setdiff(colN, tofill), perl=T))){
      v <- as.numeric( gsub('^X','', v, perl=T) )
      stopifnot(!is.na(v))
      stopifnot(v >= 0)
      stopifnot(!any(is.na(timeAvailable)))
      
      stopifnot(! v %in% timeAvailable )
      stopifnot(all(timeAvailable >= 0))
      
      timeAvailable <- sort( timeAvailable )
      if(all(timeAvailable < v)){
        x.last <- x[, paste('X', tail(timeAvailable, 1), sep='')]
        colnames(x.last) <- paste0("X", v)
        return( x.last ) # using the last time point
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
        
      }
      leftProb <- as.numeric(unlist(leftProb))
      (leftProb * (rightTime -v)  + rightProb * ( v - leftTime ) )/(rightTime - leftTime)
      
      #this is way too slow
      #return( apply(cbind(leftProb, rightProb), 1, function(p){approx(x=c(leftTime, rightTime), y=p, xout=v )$y}) )
      
      
    })
    
    #  colnames(fillin ) <- tofill
    stopifnot(all(colnames(fillin) == tofill))
    x <- cbind(x, fillin)
    
  }
  
  if(method == 'previousEstimate'){
    fillin<- sapply(tofill, function(v, timeAvailable = as.numeric( gsub('^X','', setdiff(colN, tofill), perl=T) ) ){
      v <- as.numeric( gsub('^X','', v, perl=T) )
      stopifnot(!is.na(v))
      stopifnot(!any(is.na(timeAvailable)))
      
      stopifnot(! v %in% timeAvailable )
      timeAvailable <- sort(timeAvailable)
      
      if(all(timeAvailable > v)) return('one') # earlier than any observed event time, the survival prob will be set to 1
      i <- max( which(timeAvailable < v)  )
      v1 <- paste('X', timeAvailable[i], sep='')
      browser()
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
  
  
  
  #stopifnot( all(cols %in% colnames(x)) )
  #  x[, c(cols, toaddBack),drop=F]
  x.interpolate <- x[, cols, drop=F]
  return(x.interpolate)
}

aggBrier <- function(pred, observed, outcome.var, batchVar='seed', idVar=''){
  require(pec)
  stopifnot(length(outcome.var) == 2)
  #print(paste('variable', outcome.var[1], 'from',deparse(substitute(observed)), 'will be treated as TTE'))
  #print(paste('variable', outcome.var[2], 'from',deparse(substitute(observed)), 'will be treated as event indicator (1 as event, 0 as censor)'))
  
  form_censor <- paste("Surv(",outcome.var[1],",",outcome.var[2],")~1")
  #print(paste('censoring model for IPCW', form_censor))
  form_censor <- as.formula(form_censor)
  
  stopifnot(all(outcome.var %in% colnames(observed)))
  stopifnot(all(c(batchVar ) %in% colnames(pred)))
  stopifnot(nchar(idVar) > 0)
  
  ub <- unique(pred[, batchVar])
  un <- unique(pred[, 'size'])
  
  time.c <- grep('^X[\\d\\.]+$', colnames(pred), value=T, perl=T)
  time.n <- as.numeric( gsub('X', '', time.c) )
  stopifnot(!any(is.na(time.n)))
  
  if(nchar(idVar) > 0){
    stopifnot(idVar %in% colnames(observed))
    stopifnot(idVar %in% colnames(pred))
    stopifnot(!any(duplicated(observed[, idVar])))
    row.names(observed) <- as.character(observed[, idVar])
  }  
  
  if(time.n[1] == 0){
    print('input matrix contains a column for time 0')
    stopifnot(all(pred[, time.c[1]] == 1))
  }
  
  r2 <- sapply( ub, function(b){
    r1 <- sapply(un, function(s){
      p1 <- pred[pred[, batchVar] %in% b & pred[, 'size'] %in% s,  ]
      
      if(nrow(p1) != nrow(observed)){
        print(paste('there are', nrow(p1),'predictions for', nrow(observed),'observations from batch', b,'using', s,'features'))
        return(NULL)
      }
      
      if(nchar(idVar) > 0){
        stopifnot(!any(duplicated(p1[, idVar])))
        observed <- observed[p1[, idVar], ]
        stopifnot(all( p1[, idVar] == observed[, idVar]))
      }
      
      p2 <- data.matrix(p1[, time.c[order(time.n)]] )
      if(time.n[1] != 0) {
        p2 <- cbind(1, p2)
      }
      
      brier <-  pec(object=list('prevalidationPrediction' = p2), formula=form_censor, data=observed[, outcome.var], exact=F, cens.model="marginal",splitMethod="none",B=0, verbose=F, times=time.n)
      
      #brier <- pec::pec(survPredictions, formula, data=data, times=times, exact=FALSE, verbose=FALSE) #from Harbron, Chris
      # as there are few time points with survival prediction than observed event time, specifying exact =T will lead to an error without or without providing time.n. This is different from the manual: If times are given and exact=TRUE then the times are merged with the unique values of the response variable. 
      # after providing exact =F,  not providing times=time.n will also lead to an error message.
      
      
      #specify times? need zero time point
      ibrier <-  crps(brier, times=max(time.n), start = min( time.n[time.n>0] )  )['prevalidationPrediction', 1]
      
      
    })
    
    if(all(is.null(unlist(r1)))) return(NULL)
    names(r1) <- as.character(un)
    r1
  })
  
  if(class(r2) %in% c('matrix','data.frame')) colnames(r2) <- as.character(ub)
  r2
}
