checkRetGroupSetParams <- function(params) {

  if(params$retcorMethod == "obiwarp") {
    quantitative_parameters <- 
      c("profStep", "gapInit", "gapExtend", "response", "factorDiag", 
                                 "factorGap")
    qualitative_parameters <- 
      c("distFunc", "plottype", "localAlignment", "center")
    unsupported_parameters <- 
      c("col", "ty", "initPenalty")
  } else {
    if(params$retcorMethod == "loess") {
      quantitative_parameters <- c("missing", "extra", "span")
      qualitative_parameters <- c("smooth", "family", "plottype")
      unsupported_parameters <- c("col", "ty")
    }    
  }
  
  quantitative_parameters <- 
    c(quantitative_parameters, "minfrac", "minsamp", "bw", "mzwid", "max")
  unsupported_parameters <- 
    c(unsupported_parameters, "sleep")
  
  checkParams(params, quantitative_parameters, qualitative_parameters, 
              unsupported_parameters)
}


getDefaultRetCorCenterSample <- function(xset) {
  ret <- NULL
  for(i in 1:length(filepaths(xset))) {
    ret <- c(ret, sum(peaks(xset)[,"sample"] == i))
  }
  return(which.max(ret))
}


getDefaultRetGroupStartingParams <-
  function(retcorMethod = c("obiwarp", "loess", "none"), 
           distfunc=c("cor_opt", "cor", "cov", "prd", "euc"), 
           high_resolution=TRUE) {

  retcorMethod <- match.arg(retcorMethod)
  ret <- NULL  
  
  if(retcorMethod == "obiwarp") {
    distfunc <- match.arg(distfunc)
    if(distfunc=="cor")
      ret <- (list(distFunc="cor", 
                   gapInit=c(0.0, 0.4), 
                   gapExtend=c(2.1, 2.7)))
	  if(distfunc=="cor_opt")
	    ret <- (list(distFunc="cor_opt", 
	                 gapInit=c(0.0, 0.4), 
	                 gapExtend=c(2.1, 2.7)))
    if(distfunc=="cov")
	    ret <- (list(distFunc="cov", 
	                 gapInit=c(0.0, 0.4), 
	                 gapExtend=c(11.4, 12.0)))
	  if(distfunc=="prd")
	    ret <- (list(distFunc="prd", 
	                 gapInit=c(0.0, 0.4), 
	                 gapExtend=c(7.5, 8.1)))
	  if(distfunc=="euc")
	    ret <- (list(distFunc="euc", 
	                 gapInit=c(0.7, 1.1), 
	                 gapExtend=c(1.5, 2.1)))

    ret$profStep <- c(0.7, 1)
    ret$plottype <- "none"
    ret$response <- 1
    ret$factorDiag <- 2
    ret$factorGap <- 1
    ret$localAlignment <- 0
    #ret$initPenalty <- 0
	  
  } else {
    if(retcorMethod == "loess") {
      ret <- list()
      ret$missing <- c(1,3)
      ret$extra <- c(1,3)
	    ret$span <- c(0.1, 0.3)
	    ret$smooth <- "loess"
	    ret$family <- "gaussian"
	    ret$plottype <- "none"
	  #ret$col <- NULL
	  #ret$ty <- NULL
    } else { #dont correct retention time
      ret <- list()
    }
  }
  
  ret$retcorMethod <- retcorMethod

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


getNormalizedResponse <- function(response) {

  #good_groups <- sapply(response, "[[", "good_groups")
  #bad_groups <- sapply(response, "[[", "bad_groups")
  #bad_groups_bool <- bad_groups == 0
  #bad_groups[bad_groups_bool] <- 1
  #good_groups[bad_groups_bool] <- good_groups[bad_groups_bool] + 1
  #group_ratio <- good_groups ^ 2 / bad_groups 
  GS <- sapply(response, "[[", "GS")
  RCS <- sapply(response, "[[", "RCS")
  
  #give penalty when retcor failed
  RCS_penalty <- 1/sapply(response, "[[", "retcor_done")
  RCS <- RCS/RCS_penalty
  
  #normalize
  norm_GS <- (GS - min(GS)) / (max(GS) - min(GS))  
  norm_RCS <- (RCS - min(RCS)) / (max(RCS) - min(RCS))
  norm_GS[is.na(norm_GS)] <- 0
  norm_RCS[is.na(norm_RCS)] <- 0
  
  return(norm_GS+norm_RCS)

}


getRGTVValues <- function(xset, exp_index=1, retcor_penalty=1) {

  relative_rt_diff <- c()
  
  if(nrow(xcms::groups(xset)) > 0) {
    for(i in 1:nrow(xcms::groups(xset))) {
      feature_rtmed <- xcms::groups(xset)[i, "rtmed"]
	    relative_rt_diff <- 
	      c(relative_rt_diff, 
	        mean(abs(feature_rtmed - 
	                   peaks(xset)[groupidx(xset)[[i]], "rt"]) / feature_rtmed))
    }
    good_groups <- 
      sum(unlist(lapply(X=groupidx(xset), FUN = function(x, xset) {
        ifelse(length(unique(peaks(xset)[x,"sample"])) == 
                 length(filepaths(xset)) & 
                 length(peaks(xset)[x,"sample"]) == 
                 length(filepaths(xset)), 1, 0)
      }, xset)))
    bad_groups <- nrow(xcms::groups(xset)) - good_groups
  } else {
    relative_rt_diff <- 1
    good_groups <- 0
    bad_groups <- 0   
  }
  
  tmp_good_groups <- good_groups + ifelse(bad_groups==0, 1, 0)
  tmp_bad_groups <- bad_groups + ifelse(bad_groups==0, 1, 0)
  
  ARTS <- (mean(relative_rt_diff)) * retcor_penalty
  
  ret <- list(exp_index   = exp_index, 
              good_groups = good_groups, 
              bad_groups  = bad_groups, 
              GS          = tmp_good_groups^2/tmp_bad_groups, 
              RCS         = 1/ARTS)
  
  ret$retcor_done = retcor_penalty        
  
  return(ret)  
}


optimizeRetGroup <- 
  function(xset, 
           params=getDefaultRetGroupStartingParams(), 
           nSlaves=4, 
           subdir="IPO") {
                                                 
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
    message("\n\n")
    message("starting new DoE with:\n")
    message(paste(rbind(paste(names(params), sep="", ": "), 
                        paste(params, sep="", "\n")),
                  sep=""))
        
    retcor_result <- 
      retGroupCalcExperimentsCluster(params, xset, nSlaves)  
                       
    retcor_result <- 
      retGroupExperimentStatistic(retcor_result, subdir, iterator, xset)
    
    history[[iterator]] <- retcor_result 
         
    params <- retcor_result$params$to_optimize 
	
    if(!RCSandGSIncreased(history)) {
      message("no increase stopping")
      
      history$best_settings <- history[[(length(history)-1)]]$best_settings	 
	  
      return(history)
    }
             
    for(i in 1:length(params)) {

      parameter_setting <- retcor_result$max_settings[i+1]
      bounds <- params[[i]]
      curParam <- names(params)[i]
		
	    min_bound <- ifelse(curParam == "profStep",0.3, 
                   ifelse(curParam == "mzwid",0.0001,
                   ifelse(curParam == "bw",0.25,
                   ifelse(curParam == "span",0.001,0))))		
      
      step_factor <- 
        ifelse(is.na(parameter_setting), 1.2, 
               ifelse((abs(parameter_setting) <= best_range), 0.8, 1))
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
           #if parameter are within range, decrease parameter-range
	         if(abs(parameter_setting) < best_range) { 
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
      
       if(curParam == "extra" | curParam == "missing") {
         params[[i]] <- round(params[[i]])
       }
      
       if(curParam == "profStep" | curParam == "minfrac") {
         if(params[[i]][2] > 1) {  # 1 is max value for profStep
           params[[i]] <- round(c(1-(diff(params[[i]])*0.8), 1),2)
           message(paste("profStep or minfrac greater 1, decreasing to", 
                         paste(params[[i]], collapse = " and ")))  
         }      
       }	
	  }
	
	  params <- attachList(params, retcor_result$params$no_optimization)
              
    iterator <- iterator + 1                 
  }
  return(history)
}

optimizeRetGroupSlaveCluster <- function(task, xset, parameters) {
      
  processedData <- retcorGroup(xset, parameters, task)

  result <- 
    getRGTVValues(processedData$xset, task, processedData$retcor_failed)
	return(result)

  
}



RCSandGSIncreased <- function(history) {

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
  
  if((cur_tv$good_groups^2/cur_tv$bad_groups <= 
      prev_tv$good_groups^2/prev_tv$bad_groups) | 
     (cur_tv$RCS <= prev_tv$RCS))
    return(FALSE)
    
  return(TRUE)

}


retGroupCalcExperimentsCluster <- function(params, xset, nSlaves=4) {

  typ_params <- typeCastParams(params)
  
  if(length(typ_params$to_optimize)>1) {
    design <- getCcdParameter(typ_params$to_optimize)    
    parameters <- decode.data(design) 
  } else {
    design <- data.frame(run.order=1:9, a=seq(-1,1,0.25))
    colnames(design)[2] <- names(typ_params$to_optimize)
    parameters <- design
    parameters[,2] <- 
      seq(min(typ_params$to_optimize[[1]]), 
          max(typ_params$to_optimize[[1]]), 
          diff(typ_params$to_optimize[[1]])/8)
  }  
  
  tasks <- as.list(1:nrow(design))      
  parameters <- combineParams(parameters, typ_params$no_optimization)
  
  if(nSlaves > 1) {
    # unload snow (if loaded) to prevent conflicts with usage of 
    # package 'parallel'
    if('snow' %in% rownames(installed.packages()))
      unloadNamespace("snow")
    
    cl_type<-getClusterType()
    cl <- parallel::makeCluster(nSlaves, type = cl_type) #, outfile="log.txt")
    #exporting all functions to cluster but only calcRGTV is needed
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    if(identical(cl_type,"PSOCK")) {
      message("Using PSOCK type cluster, this increases memory requirements.")
      message("Reduce number of slaves if you have out of memory errors.")
      message("Exporting variables to cluster...")
      parallel::clusterExport(cl, ex)
    }
    result <- parallel::parLapply(cl, tasks, optimizeRetGroupSlaveCluster, xset, 
                                  parameters)
    parallel::stopCluster(cl)
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

retGroupExperimentStatistic <- function(retcor_result, subdir, iterator, xset) {

  params <- retcor_result$params
  resp <- getNormalizedResponse(retcor_result$response)
  
  model <- createModel(retcor_result$design, params$to_optimize, resp)
  retcor_result$model <- model                  
  
  max_settings <- getMaximumLevels(retcor_result$model)
  tmp <- max_settings[1,-1]
  tmp[is.na(tmp)] <- 1
  if(!is.null(subdir) & length(tmp) > 1)
    plotContours(retcor_result$model, 
                 tmp, 
                 paste(subdir, "/retgroup_rsm_", iterator, sep = ""))
      
  parameters <- as.list(decodeAll(max_settings[-1], params$to_optimize)) 
  parameters <- combineParams(parameters, params$no_optimization)
  xset_tmp <- xset
  
  exp_index <- 1  
  processedData <- retcorGroup(xset_tmp, parameters, exp_index)

  tv <- 
    getRGTVValues(processedData$xset, exp_index, processedData$retcor_failed)

  retcor_result$max_settings <- max_settings
  retcor_result$target_value <- tv   
  retcor_result$xset <- processedData$xset_tmp
  retcor_result$best_settings <- parameters

  return(retcor_result)

}

retcorGroup <- function(xset, parameters, exp_index=1) {
  do_retcor <- parameters$retcorMethod[exp_index] != "none"
  
  retcor_failed = ifelse(do_retcor, 1.1, 1) 
  
  if(parameters$retcorMethod[exp_index] == "loess") {
    try(
      xset <- group(
        xset, 
        method  = "density", 
        bw      = parameters$bw[exp_index], 
        mzwid   = parameters$mzwid[exp_index], 
        minfrac = parameters$minfrac[exp_index], 
        minsamp = parameters$minsamp[exp_index], 
        max     = parameters$max[exp_index])
      )
    
    try(
      retcor_failed <- retcor(
        xset, 
        method   = "loess", 
        plottype = parameters$plottype[exp_index], 
        family   = parameters$family[exp_index],
        missing  = parameters$missing[exp_index], 
        extra    = parameters$extra[exp_index], 
        span     = parameters$span[exp_index])
      )  		  
    
    if(!is.numeric(retcor_failed)) {
      xset <- retcor_failed
      retcor_failed=1
    } 
  }
  
  if(parameters$retcorMethod[exp_index] == "obiwarp") {
    try(
      retcor_failed <- 
        retcor(xset, 
               method         = "obiwarp", 
               plottype       = parameters$plottype[exp_index], 
               distFunc       = parameters$distFunc[exp_index],
               profStep       = parameters$profStep[exp_index], 
               center         = parameters$center[exp_index], 
               response       = parameters$response[exp_index], 
               gapInit        = parameters$gapInit[exp_index], 
               gapExtend      = parameters$gapExtend[exp_index],
               factorDiag     = parameters$factorDiag[exp_index],
               factorGap      = parameters$factorGap[exp_index], 
               localAlignment = parameters$localAlignment[exp_index])
      )  	
    
    if(!is.numeric(retcor_failed)) {
      xset <- retcor_failed
      retcor_failed = 1
    }       
  } 
  
  try(
    xset <- group(
      xset, 
      method  = "density", 
      bw      = parameters$bw[exp_index], 
      mzwid   = parameters$mzwid[exp_index], 
      minfrac = parameters$minfrac[exp_index], 
      minsamp = parameters$minsamp[exp_index], 
      max     = parameters$max[exp_index])
    )
  
  
  return(list(xset = xset, retcor_failed = retcor_failed))
}


