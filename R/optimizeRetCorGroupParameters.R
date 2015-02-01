checkRetGroupSetParams <-
function(params) {

  quantitative_parameters <- c("profStep", "gapInit", "gapExtend", "response", "factorDiag", "factorGap", "minfrac", "minsamp", "bw", "mzwid", "max")
  qualitative_parameters <- c("distFunc", "plottype", "localAlignment", "center")
  unsupported_parameters <- c("col", "ty", "initPenalty", "sleep")

  checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
}


getDefaultRetCorCenterSample <-
function(xset) {
  ret <- NULL
  for(i in 1:length(xset@filepaths)) {
    ret <- c(ret, sum(xset@peaks[,"sample"] == i))
  }
  return(which.max(ret))
}


getDefaultRetGroupStartingParams <-
function(distfunc="cor_opt", high_resolution=TRUE) {

  ret <- NULL
  if(!is.null(distfunc)) {
    if(distfunc=="cor")
      ret <- (list(distFunc="cor", gapInit=c(0.0, 0.4), gapExtend=c(2.1, 2.7)))
	if(distfunc=="cor_opt")
	  ret <- (list(distFunc="cor_opt", gapInit=c(0.0, 0.4), gapExtend=c(2.1, 2.7)))
    if(distfunc=="cov")
	  ret <- (list(distFunc="cov", gapInit=c(0.0, 0.4), gapExtend=c(11.4, 12.0)))
	if(distfunc=="prd")
	  ret <- (list(distFunc="prd", gapInit=c(0.0, 0.4), gapExtend=c(7.5, 8.1)))
	if(distfunc=="euc")
	  ret <- (list(distFunc="euc", gapInit=c(0.7, 1.1), gapExtend=c(1.5, 2.1)))

    ret$profStep <- c(0.7, 1)
    ret$plottype <- "none"
    ret$response <- 1
    ret$factorDiag <- 2
    ret$factorGap <- 1
    ret$localAlignment <- 0
    #ret$initPenalty <- 0
	  
  } else {
	ret <- list()
  }

	#grouping parameter
  ret$bw <- c(22,38)
  ret$minfrac <- c(0.3, 0.7)
  ret$mzwid <- c(0.15, 0.35)
  if(high_resolution)
    ret$mzwid <- c(0.015, 0.035)
  ret$minsamp <- 1
  ret$max <- 50
  #ret$sleep <- 0
	  
  return(ret)

}


getNormalizedResponse <-
function(response) {

  good_groups <- sapply(response, "[[", "good_groups")
  bad_groups <- sapply(response, "[[", "bad_groups")
  bad_groups_bool <- bad_groups == 0
  bad_groups[bad_groups_bool] <- 1
  good_groups[bad_groups_bool] <- good_groups[bad_groups_bool] + 1
  group_ratio <- good_groups ^ 2 / bad_groups 
  ARTS <- 1/sapply(response, "[[", "mean_rel_rt_diff")
  
  #give penalty when retcor failed
  ARTS_penalty <- 1/sapply(response, "[[", "retcor_done")
  ARTS <- ARTS/ARTS_penalty
  
  #normalize
  norm_group_ratio <- (group_ratio - min(group_ratio)) / (max(group_ratio) - min(group_ratio))
  norm_ARTS <- (ARTS - min(ARTS)) / (max(ARTS) - min(ARTS))
  
  return(norm_group_ratio+norm_ARTS)

}


getRGTVValues <-
function(xset, exp_index=1, retcor_penalty=1) {
  features <- nrow(xset@groups)

  good_groups <- sum(unlist(lapply(X=xset@groupidx, FUN=function(x, xset) {ifelse(length(unique(xset@peaks[x,"sample"]))==length(xset@filepaths) & 
                                   length(xset@peaks[x,"sample"])==length(xset@filepaths),1,0)}, xset)))
  bad_groups <- nrow(xset@groups) - good_groups

  relative_rt_diff <- c()
  
  if(nrow(xset@groups) > 0) {
    for(i in 1:nrow(xset@groups)) {
      feature_rtmed <- xset@groups[i, "rtmed"]
	    relative_rt_diff <- c(relative_rt_diff, mean(abs(feature_rtmed - xset@peaks[xset@groupidx[[i]], "rt"])/feature_rtmed))
    }
  } else {
    relative_rt_diff <- 1
  }
  
  ARTS <- (mean(relative_rt_diff)) * retcor_penalty
  
  ret <- list(exp_index=exp_index, good_groups=good_groups, bad_groups=bad_groups, mean_rel_rt_diff=ARTS)
  
  ret$retcor_done = retcor_penalty        
  
  return(ret)  
}


optimizeRetGroup <-
function(xset, params=getDefaultRetGroupStartingParams(), nSlaves=4, subdir="IPO") {
                                                 
  iterator = 1 
  history <- list()  
  best_range <- 0.25

  if(!is.null(subdir))
    if(!file.exists(subdir))
      dir.create(subdir)
	
  if(is.null(params$center))
    params$center <- getDefaultRetCorCenterSample(xset) 

  checkRetGroupSetParams(params)	
    
  while(iterator < 50) {
    cat("\n")
    cat("\n")
    cat("\n")
    cat("starting new DoE with:\n")
    print(params)
        
    retcor_result <- retGroupCalcExperimentsCluster(params, xset, nSlaves)  
                       
    retcor_result <- retGroupExperimentStatistic(retcor_result, subdir, iterator, xset)
    
    history[[iterator]] <- retcor_result 
         
    params <- retcor_result$params$to_optimize 
	#all settings must be within range, not like xcmsSet-Optimization
    #found_good_settings <- foundGoodSettings(rep(TRUE, length(retcor_result$parameter_significance)), retcor_result$eigenvalues, retcor_result$max_settings[-1])
    
    if(!RCSandGSIncreased(history)) {
      print("no increase stopping")
      
      history$best_settings <- history[[(length(history)-1)]]$best_settings	 
	  
      return(history)       

    }
             
    for(i in 1:length(params)) {

      parameter_setting <- retcor_result$max_settings[i+1]
      bounds <- params[[i]] 
		
	    min_bound <- ifelse(names(params)[i] == "profStep",0.3, ifelse(names(params)[i] == "mzwid",0.0001,ifelse(names(params)[i] == "bw",0.25,0)))		
      
      step_factor <- ifelse(is.na(parameter_setting), 1.2, ifelse((abs(parameter_setting) < best_range), 0.8, 1))
      step <- (diff(bounds) / 2) * step_factor
      
      if(is.na(parameter_setting)) {      
		    new_center <- decode(0, bounds)
		    step <- (diff(bounds) / 2) * 1.2
		    if(new_center > step) {
		      new_bounds <- c(new_center - step, new_center + step) 
		    } else {
		      new_bounds <- c(min_bound, 2*step + min_bound)      
		    }

		    names(new_bounds) <- NULL
            params[[i]] <- new_bounds
			
		  } else {		  
		    if(bounds[1] > min_bound | parameter_setting > -1) {
          #center around optimum
          new_center <- decode(parameter_setting, bounds)
          step <- diff(bounds) / 2
	      if(abs(parameter_setting) < best_range) { #if parameter are within range, decrease parameter-range
	        step <- step * 0.8
        }

			  if(new_center > step) {
          new_bounds <- c(new_center - step, new_center + step) 
			  } else {
		      new_bounds <- c(min_bound, 2*step + min_bound)  
        }			  

		    if(parameter_setting == -1 || parameter_setting == 1) {
		      step <- step / 5
		      new_bounds[1] <- max(min_bound, new_bounds[1] - step)
		      new_bounds[2] <- new_bounds[2] + step
		    } 
			#print(paste("new_bounds", paste(new_bounds, collapse=" ")))
            names(new_bounds) <- NULL
			
            params[[i]] <- new_bounds
          } else {
		        step <- diff(bounds) / 2
            params[[i]] <- c(min_bound, 2*step + min_bound)
		      }		  
      }
      
      if(names(params)[i] == "profStep" | names(params)[i] == "minfrac") {
        if(params[[i]][2] > 1) {  # 1 is max value for profStep
          params[[i]] <- round(c(1-(diff(params[[i]])*0.8), 1),2)
          print(paste("profStep or minfrac greater 1, decreasing to", paste(params[[i]], collapse = " and ")))  
        }      
      }
		
	}
	
	params <- attachList(params, retcor_result$params$no_optimization)
              
    iterator <- iterator + 1
                 
  }
  return(history)

}

optimizeRetGroupSlaveCluster <-
function(task, xset, parameters) {
	#print(parameters)
    #library(xcms)
    #if (tag == 1) {
      exp_index <- task

      do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && is.null(parameters$gapInit) && is.null(parameters$gapExtend)
	                && is.null(parameters$plottype) && is.null(parameters$response) 
					&& is.null(parameters$factorDiag) && is.null(parameters$factorGap) && is.null(parameters$localAlignment) 
					&& is.null(parameters$initPenalty)) #&& is.null(parameters$center)) 
					
      retcor_failed = ifelse(do_retcor, 1.1, 1)  
  
      if(do_retcor) {
        try(retcor_failed <- retcor(xset, method="obiwarp", plottype=parameters$plottype[exp_index], distFunc=parameters$distFunc[exp_index],
                             profStep=parameters$profStep[exp_index], center=parameters$center[exp_index], response=parameters$response[exp_index], 
							 gapInit=parameters$gapInit[exp_index], gapExtend=parameters$gapExtend[exp_index],
							 factorDiag=parameters$factorDiag[exp_index], factorGap=parameters$factorGap[exp_index], 
							 localAlignment=parameters$localAlignment[exp_index]))
  	
	
        if(!is.numeric(retcor_failed)) {
          xset <- retcor_failed
          retcor_failed=1
        } 
      }
	  
      minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])      
      try(xset <- group(xset, method="density", bw=parameters$bw[exp_index], mzwid=parameters$mzwid[exp_index], minfrac=minfrac, 
                  minsamp=parameters$minsamp[exp_index], max=parameters$max[exp_index]))	 
	  #print(xset)
      result <- getRGTVValues(xset, exp_index, retcor_failed)
	  return(result)
      #mpi.send.Robj(result,0,2)

    #} else if (tag == 2) {
    #  done <- 1
    #}
    # Else ignore the message or report an error
  #}

  # Tell master that this slave is exiting.  Send master an exiting message
  #mpi.send.Robj(junk,0,3) 
  
}


# optimizeRetGroupSlave <-
# function() {
  # junk <- 0
  # done <- 0
  
  # while (done != 1) {
  #  #Signal being ready to receive a new task
    # mpi.send.Robj(junk,0,1)

  #  #Receive a task
    # task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    # task_info <- mpi.get.sourcetag()
    # tag <- task_info[2]
    # print(paste("tag", tag, "task", task))
	
	# print(parameters)
    
    # if (tag == 1) {
      # exp_index <- task

      # do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && is.null(parameters$gapInit) && is.null(parameters$gapExtend)
	                # && is.null(parameters$plottype) && is.null(parameters$response) 
					# && is.null(parameters$factorDiag) && is.null(parameters$factorGap) && is.null(parameters$localAlignment) 
					# && is.null(parameters$initPenalty)) #&& is.null(parameters$center)) 
					
      # retcor_failed = ifelse(do_retcor, 1.1, 1)  
  
      # if(do_retcor) {
        # try(retcor_failed <- retcor(xset, method="obiwarp", plottype=parameters$plottype[exp_index], distFunc=parameters$distFunc[exp_index],
                             # profStep=parameters$profStep[exp_index], center=parameters$center[exp_index], response=parameters$response[exp_index], 
							 # gapInit=parameters$gapInit[exp_index], gapExtend=parameters$gapExtend[exp_index],
							 # factorDiag=parameters$factorDiag[exp_index], factorGap=parameters$factorGap[exp_index], 
							 # localAlignment=parameters$localAlignment[exp_index]))
  	
	
        # if(!is.numeric(retcor_failed)) {
          # xset <- retcor_failed
          # retcor_failed=1
        # } 
      # }
	  
      # minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])      
      # try(xset <- group(xset, method="density", bw=parameters$bw[exp_index], mzwid=parameters$mzwid[exp_index], minfrac=minfrac, 
                  # minsamp=parameters$minsamp[exp_index], max=parameters$max[exp_index]))	 
 
      # result <- getRGTVValues(xset, exp_index, retcor_failed)

      # mpi.send.Robj(result,0,2)

    # } else if (tag == 2) {
      # done <- 1
    # }
  #  #Else ignore the message or report an error
  # }

#  #Tell master that this slave is exiting.  Send master an exiting message
  # mpi.send.Robj(junk,0,3) 
  
# }


RCSandGSIncreased <-
function(history) {

  index = length(history)
  if(index < 2)
    return(TRUE)
    
  prev_tv <- history[[index-1]]$target_value
  cur_tv <- history[[index]]$target_value
  
  if(cur_tv$bad_groups == 0) {
    cur_tv$bad_groups = 1
    cur_tv$good_groups = cur_tv$good_groups + 1
  }
  
  if(prev_tv$bad_groups == 0) {
    prev_tv$bad_groups = 1
    prev_tv$good_groups = prev_tv$good_groups + 1
  }
  
  if((cur_tv$good_groups^2/cur_tv$bad_groups <= prev_tv$good_groups^2/prev_tv$bad_groups) | (cur_tv$mean_rel_rt_diff >= prev_tv$mean_rel_rt_diff))
    return(FALSE)
    
  return(TRUE)

}


retGroupCalcExperimentsCluster <-
function(params, xset, nSlaves=4) {

  typ_params <- typeCastParams(params)
  
  #if(length(typ_params$to_optimize) > 2) {
  #  design <- getBbdParameter(typ_params$to_optimize) 
  #} else {
    design <- getCcdParameter(typ_params$to_optimize) 
  #}	
 
  parameters <- decode.data(design)	
  tasks <- as.list(1:nrow(design))      
  parameters <- combineParams(parameters, typ_params$no_optimization)
  
  if(nSlaves > 1) {
    cl <- makeCluster(nSlaves, type = "PSOCK")#, outfile="log.txt")
  #exporting all functions to cluster but only calcRGTV is needed
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    clusterExport(cl, ex)
    result <- parLapply(cl, tasks, optimizeRetGroupSlaveCluster, xset, parameters)
    stopCluster(cl)
  } else {
	result <- lapply(tasks, optimizeRetGroupSlaveCluster, xset, parameters)
  }  
  #print(result) 
  response <- list()   
  for(i in 1:length(result))
    response[[result[[i]]$exp_index]] <- result[[i]]
  
  ret <- list()
  ret$params <- typ_params
  ret$design <- design
  ret$response <- response
  
  return(ret)

}

#retGroupCalcExperiments <-
#function(params, xset, nSlaves=4) {
#    
#  junk <- 0
#  closed_slaves <- 0
#  #nSlaves <- min(mpi.comm.size()-1, nSlaves)  
#  
#  typ_params <- typeCastParams(params)
#  
#  if(length(typ_params$to_optimize) > 2) {
#    design <- getBbdParameter(typ_params$to_optimize) 
#  } else {
#    design <- getCcdParameter(typ_params$to_optimize) 
#  }	
# 
#  parameters <- decode.data(design)	
#  tasks <- as.list(1:nrow(design))    
#  startSlaves(nSlaves)
#  
#  parameters <- combineParams(parameters, typ_params$no_optimization)
#  
#  sendRetGroupSlaveFunctions(parameters, xset) 
#  
#  response <- list()  
#  finished <- 0
#  while(closed_slaves < nSlaves) {
#    # Receive a message from a slave
#    message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
#    message_info <- mpi.get.sourcetag()
#    slave_id <- message_info[1]
#    tag <- message_info[2]
#      
#    if(tag == 1) {
#      if(length(tasks) > 0) {
#        mpi.send.Robj(tasks[[1]], slave_id, 1);
#        tasks[[1]] <- NULL
#      } else {
#        mpi.send.Robj(junk, slave_id, 2)
#      }
#    } else if (tag == 2) {
#	  #print(message)
#      response[[message$exp_index]] <- message
#      finished <- finished + 1     
#    } else if (tag == 3) {
#      # A slave has closed down. 
#      closed_slaves <- closed_slaves + 1
#    }
#    cat(paste("finished ", finished, " of ", nrow(design), " tasks\r", sep="")) 
#    flush.console()
#  }
#  cat("\n\r")  
#  print("done")
#
#  mpi.close.Rslaves()
#  
#  ret <- list()
#  ret$params <- typ_params
#  ret$design <- design
#  #ret$model <- model
#  ret$response <- response
#  
#  return(ret)
#
#}

retGroupExperimentStatistic <-
function(retcor_result, subdir, iterator, xset) {

  params <- retcor_result$params
  resp <- getNormalizedResponse(retcor_result$response)
  model <- createModel(retcor_result$design, params$to_optimize, resp)
  
  retcor_result$model <- model                  
  max_settings <- getMaximumExperiment(retcor_result$model)
  tmp <- max_settings[-1]
  tmp[is.na(tmp)] <- 1

  if(!is.null(subdir))
    plotContours(retcor_result$model, tmp, paste(subdir, "/retgroup_rsm_", iterator, sep=""))  
    
  parameters <- as.list(decodeAll(max_settings[-1], params$to_optimize)) 
  parameters <- combineParams(parameters, params$no_optimization)
  xset_tmp <- xset
  
  exp_index <- 1				   
  do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && is.null(parameters$gapInit) && is.null(parameters$gapExtend)
	             && is.null(parameters$plottype) && is.null(parameters$col) && is.null(parameters$ty) && is.null(parameters$response) 
				 && is.null(parameters$factorDiag) && is.null(parameters$factorGap) && is.null(parameters$localAlignment) 
				 && is.null(parameters$initPenalty)) #&& is.null(parameters$center)) 
					
      retcor_failed = ifelse(do_retcor, 1.1, 1)  
  
  if(do_retcor) {
    try(retcor_failed <- retcor(xset_tmp, method="obiwarp", plottype=parameters$plottype[exp_index], distFunc=parameters$distFunc[exp_index],
                             profStep=parameters$profStep[exp_index], center=parameters$center[exp_index],
							 response=parameters$response[exp_index], gapInit=parameters$gapInit[exp_index], 
				 	         gapExtend=parameters$gapExtend[exp_index], factorDiag=parameters$factorDiag[exp_index], 
					         factorGap=parameters$factorGap[exp_index], localAlignment=parameters$localAlignment[exp_index]))
  	
	
    if(!is.numeric(retcor_failed)) {
      xset_tmp <- retcor_failed
      retcor_failed=1
    } 
  }
	  
  minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])      
  try(xset_tmp <- group(xset_tmp, method="density", bw=parameters$bw[exp_index], mzwid=parameters$mzwid[exp_index], minfrac=minfrac, 
                  minsamp=parameters$minsamp[exp_index], max=parameters$max[exp_index]))	 
 
  tv <- getRGTVValues(xset_tmp, exp_index, retcor_failed)

  retcor_result$max_settings <- max_settings
  retcor_result$target_value <- tv   
  retcor_result$xset <- xset_tmp
  retcor_result$best_settings <- parameters

  return(retcor_result)

}


# sendRetGroupSlaveFunctions <-
# function(parameters, xset) {
  # mpi.bcast.Robj2slave(getRGTVValues)
  # mpi.bcast.cmd(slave <- mpi.comm.rank())
  # mpi.bcast.Robj2slave(parameters)
  # mpi.bcast.Robj2slave(xset)  
  # mpi.bcast.Robj2slave(optimizeRetGroupSlave)  
  # mpi.bcast.cmd("library(xcms)")
  # mpi.bcast.cmd(optimizeRetGroupSlave())
# }

