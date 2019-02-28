#' generate command strings
#' @param runSpec list, user inputs
#' @param curr_seed int, seed number to be used
#' @param cv int, cv number to be used
#' @param file_run, log filename for printing temp results while job is running
#' @param file_cmd, filename where the sbatch command is saved into
#' 
generate_cmd <- function(runSpec, curr_seed, cv, file_run, file_cmd){
  log_path <- runSpec["log_path"]

  if(grepl("(py)$", runSpec$engine)){
    scripts <- paste0("#!/bin/bash
                        #BSUB -J ",  runSpec[['project_name']],"
                        #BSUB -oo ", runSpec[['project_name']],"_%J        
                        #BSUB -eo ", runSpec[['project_name']],"_%J",
                      "\npython3 ", runSpec["engineFile"], 
                      " ", curr_seed,
                      " ", paste0(log_path, '/metainfo.txt'),
                      " ", cv)
    cat(scripts, file=paste0(log_path, "/", file_cmd))
    command <- paste0('bsub -q ', runSpec["queue_priority"], 
                      ' -n 12 -R "span[hosts=1] rusage[mem=4]" -oo ', log_path, "/", file_run, 
                      ' < ', log_path, "/", file_cmd)
    return(command)
    
  } else {
    scripts <- paste0("#!/bin/bash
                        #BSUB -J ",  runSpec[['project_name']],"
                        #BSUB -oo ", runSpec[['project_name']],"_%J        
                        #BSUB -eo ", runSpec[['project_name']],"_%J",
                      "\nR CMD BATCH   --no-save --no-restore '--args ", curr_seed, ## arg 1
                      " ", paste0(log_path, '/metainfo.Rdata'), ## arg 2
                      " ", cv, ## arg 3, CV number
                      "' ", runSpec["engineFile"], ## engine file
                      " ", log_path, "/", file_run)
    cat(scripts, file=paste0(log_path, "/", file_cmd))
    command <- paste0("bsub -q ", runSpec["queue_priority"],
                      ' -n 12 -R "span[hosts=1] rusage[mem=4]" < ', log_path, "/", file_cmd)
    
    return(command)
  }
}

#' functions for submitting jobs
#' @param runSpec list, user inputs

sbatch_submit <- function(runSpec){
  
  engine_path <- "~/EVE/eve/models/"
  #engineFile <- paste0(engine_path, runSpec$engine)
  runSpec["engineFile"] <- paste0(engine_path, runSpec$engine)
  
  log_path <- paste0(runSpec$project_home, "/log/", runSpec$project_name)
  log_path_cluster <- paste0(log_path, "/cluster_out")
  runSpec["log_path"] <- log_path
  runSpec["log_path_cluster"] <- log_path_cluster
  
  runSpec$outP <- paste0(runSpec$project_home, "/results/", runSpec$project_name)
  
  ## get current EVE's commit hash
  chash <- git2r::revparse_single(git2r::repository('~/EVE'),"HEAD")
  runSpec["CommitHash"] <- paste(chash, collapse = ",")
  
  ## if log directly not exist, create one
  if(!dir.exists(log_path)){
    dir.create(log_path, recursive=TRUE)
  }
  ## if log/cluster_out not exist, create one
  if(!dir.exists(log_path_cluster)){
    dir.create(log_path_cluster, recursive=TRUE)
  }
  setwd(log_path_cluster)
  
  ####################################
  ## save runSpec to the log folder ##
  ####################################
  ## check if LOOCV
  if(runSpec$num_CV == -1){
    runSpec$split_CVs <- TRUE
    runSpec$num_seeds <- 1
    
    ## num_CV becomes the total number of samples
    df <- read.csv(paste(runSpec$project_home, runSpec$training_data, sep='/'), 
                   as.is=T)
    runSpec$num_CV <- nrow(df)
  }
  
  # if python engine, save runSpec as txt
  if(grepl("(py)$", runSpec$engine)){
    capture.output(print(runSpec), 
                   file = paste0(log_path, '/metainfo.txt'))
    save(runSpec, file=paste0(log_path, '/metainfo.Rdata'))
    
  } else if(grepl("r$|R$", runSpec$engine)){
    save(runSpec, file=paste0(log_path, '/metainfo.Rdata'))
    
  } else {
    stop("engine is neither .r, .R or .py")
  }
  #################
  ## submit jobs ##
  #################
  
  if(runSpec$split_CVs){
    for (seed in 1:runSpec$num_seeds){
      for (cv in 1:runSpec$num_CV){
        
        curr_seed <- seed + runSpec$seed_base
        
        file_run <- paste0("run_cv", cv, "_", curr_seed, ".log")
        file_cmd <- paste0("cmd_cv", cv, "_", curr_seed, ".log")
        
        ## create batch file command
        command <- generate_cmd(runSpec, curr_seed, cv, file_run, file_cmd)
        system(paste0("module load apps/python3
                    ", command))
      }
    }
    print(paste0(runSpec$num_seeds*cv, " jobs submitted!"))
  } else {
    for (seed in 1:runSpec$num_seeds){
      
      cv <- 0 ## indicate running all CVs in one job
      curr_seed <- seed + runSpec$seed_base
      
      file_run <- paste0("run_cv", cv, "_", curr_seed, ".log")
      file_cmd <- paste0("cmd_cv", cv, "_", curr_seed, ".log")
      
      ## create batch file command
      command <- generate_cmd(runSpec, curr_seed, cv, file_run, file_cmd)
      system(paste0("module load apps/python3
                    ", command))
    }
    print(paste0(runSpec$num_seeds, " jobs submitted!"))
  }
}


# impute missing values in the full data
imputeWithSummaryStat <- function(dsin, FUN=median, flag.var = '_F'){
  stopifnot(is.data.frame(dsin))
  m1 <- sapply(names(dsin), function(x) sum(is.na(dsin[[x]])))
  m1 <- m1[m1> 0]

  print(paste('imputing', length(m1), 'variables from', deparse(substitute(dsin)), 'with per variable', deparse(substitute(FUN))))

  for (v in names(m1)) {
      stopifnot(v %in% colnames(dsin))
      sel <- is.na(dsin[[v]])
      
      if(is.numeric(dsin[[v]])){
          s1 <- FUN(dsin[!sel, v])
          if(is.integer(dsin[[v]])) {
              s1 <- as.integer(round(s1))
          }
          stopifnot(!is.na(s1))
          print(paste('imputing', sum(sel),'missing', v,'with', round( s1, 2),', estimated from', sum(!sel),'observed values'))
          dsin[sel, v] <- s1
      }else{
          if(is.factor(dsin[[v]])) dsin[[v]] <- as.character(dsin[[v]])
          if(is.character(dsin[[v]])) {
          print(paste('set', sum(sel),'missing', v,'as UNK'))
          dsin[sel, v] <- 'UNK'
          }else{
              stop(paste('do not know how to handle variable', v,'as a', class(dsin[[v]])))
          }
      }
      
      if(!is.na(flag.var) && nchar(flag.var) > 0 ){ # create a flag variable
          f1 <- paste(v, flag.var, sep='')
          assertthat::assert_that(!f1 %in% colnames(dsin), msg = paste(f1, 'a flag variable to indicate missingness is already in the data matrix'))
          
          dsin[[f1]] <- as.integer(sel)
      }
      
  }
  dsin
}