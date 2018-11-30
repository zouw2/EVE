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
                      "\npython ", runSpec["engineFile"], 
                      " ", curr_seed,
                      " ", paste0(log_path, '/metainfo.txt'),
                      " ", cv)
    cat(scripts, file=paste0(log_path, "/", file_cmd))
    command <- paste0('bsub -q ', runSpec["queue_priority"], 
                      ' -n 6 -R "span[hosts=1] rusage[mem=12]" -oo ', log_path, "/", file_run, 
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
    command <- paste0("bsub -q ", runSpec["queue_priority"], " < ", log_path, "/", file_cmd)
    
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
    system("module load apps/python")
    
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
        system(command)
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
      system(command)
    }
    print(paste0(runSpec$num_seeds, " jobs submitted!"))
  }
}
