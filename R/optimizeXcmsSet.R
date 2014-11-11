optimizeXcmsSet <-
function(params=getDefaultStartingXcmsParams(), n_slaves=4, subdir="IPO") { #ppm=5, rt_diff=0.02, n_slaves=4, subdir="IPO") {

  checkXcmsSetParams(params)
   
  example_sample <- list.files(full.names=TRUE, ignore.case=TRUE, recursive=TRUE, pattern="(*.mzX?ML$)") 
  if(length(example_sample)==0)
    stop("no files in directory, stopping!")
	
  centWave <- is.null(params$fwhm)  
  
  history <- list()
  iterator = 1 
  best_range <- 0.25

  if(!file.exists(subdir))
    dir.create(file.path(getwd(), subdir))
    
  while(iterator < 50) {
    cat("\n")
    cat("\n")
    cat("\n")
    cat("starting new DoE with:\n")
    print(params)
        
    xcms_result <- xcmsSetExperiments(example_sample, params, n_slaves) 
#                       ppm, rt_diff, n_slaves)                        
                       
    xcms_result <- xcmsSetStatistic(xcms_result, subdir, iterator)
    history[[iterator]] <- xcms_result     
    params <- xcms_result$params 
    
    if(!resultIncreased(history)) {
      cat("no increase, stopping")
      maxima <- 0
	    max_index <- 1
	    for(i in 1:length(history)) {
	      if(history[[i]]$max_settings[1] > maxima) {
		      maxima <- history[[i]]$max_settings[1]
		      max_index <- i
		    }
	    }
   
        xcms_parameters <- as.list(decodeAll(history[[max_index]]$max_settings[-1], history[[max_index]]$params$to_optimize))      
	    xcms_parameters <- combineParams(xcms_parameters, params$no_optimization)
      
	    if(!is.list(xcms_parameters))
	      xcms_parameters <- as.list(xcms_parameters)
		
	    best_settings <- list()
      best_settings$parameters <- xcms_parameters
      library(xcms)
	  library(Rmpi)
      
      
      if(centWave) {		
        xset <- xcmsSet(files=example_sample, method="centWave", 
                  peakwidth=c(xcms_parameters$min_peakwidth, xcms_parameters$max_peakwidth),
                  ppm=xcms_parameters$ppm, noise=xcms_parameters$noise, 
				          snthresh=xcms_parameters$snthresh, mzdiff=xcms_parameters$mzdiff,
				          prefilter=c(xcms_parameters$prefilter, xcms_parameters$value_of_prefilter),
				          mzCenterFun=xcms_parameters$mzCenterFun, integrate=xcms_parameters$integrate,
				          fitgauss=xcms_parameters$fitgauss, verbose.columns=xcms_parameters$verbose.columns, 
                  nSlaves=n_slaves)
      } else {
        xset <- xcmsSet(files=example_sample, method="matchedFilter", 
                  fwhm=xcms_parameters$fwhm, snthresh=xcms_parameters$snthresh,
                  step=xcms_parameters$step, steps=xcms_parameters$steps,
                  sigma=xcms_parameters$sigma, max=xcms_parameters$max, 
                  mzdiff=xcms_parameters$mzdiff, index=xcms_parameters$index) 
                
                
      }
                
	  best_settings$xset <- xset
      target_value <- calcPPS(xset)#, ppm, rt_diff)
      best_settings$result <- target_value
      history$best_settings <- best_settings
      
      print(best_settings)
      
      return(history)   
    }

         
    for(i in 1:length(params$to_optimize)) {
      parameter_setting <- xcms_result$max_settings[i+1]
      bounds <- params$to_optimize[[i]] 
      fact <- names(params$to_optimize)[i]
	  min_factor <- ifelse(fact=="min_peakwidth", 5, ifelse(fact=="mzdiff", ifelse(centWave,-100000000, 0.001), ifelse(fact=="step",0.0005,1)))

      #if the parameter is NA, we increase the range by 20%, if it was within the inner 25% of the previous range or at the minimum value we decrease the range by 20%
      step_factor <- ifelse(is.na(parameter_setting), 1.2, ifelse((abs(parameter_setting) < best_range), 0.8, ifelse(parameter_setting==-1 & decode(-1, params$to_optimize[[i]]) == min_factor,0.8,1)))
      step <- (diff(bounds) / 2) * step_factor
      
      if(is.na(parameter_setting))
        parameter_setting <- 0
      new_center <- decode(parameter_setting, bounds)
      
      if((new_center-min_factor) > step) {
        new_bounds <- c(new_center - step, new_center + step) 
      } else {
        new_bounds <- c(min_factor, 2*step+min_factor) 
      }      
		  
	  names(new_bounds) <- NULL         
          
      if(names(params$to_optimize)[i] != "mzdiff" & names(params$to_optimize)[i] != "step")
        params$to_optimize[[i]] <- round(new_bounds, 0)
      else 
        params$to_optimize[[i]] <- new_bounds
    } 
    
    if(centWave) {
    #checking peakwidths plausiability
	  if(!is.null(params$to_optimize$min_peakwidth) | !is.null(params$to_optimize$max_peakwidth)) {
	    pw_min <- ifelse(is.null(params$to_optimize$min_peakwidth), params$no_optimization$min_peakwidth, max(params$to_optimize$min_peakwidth))
		  pw_max <- ifelse(is.null(params$to_optimize$max_peakwidth), params$no_optimization$max_peakwidth, min(params$to_optimize$max_peakwidth))
        if(pw_min >= pw_max) {
          additional <- abs(pw_min-pw_max) + 1
          if(!is.null(params$to_optimize$max_peakwidth)) {		  
            params$to_optimize$max_peakwidth <- params$to_optimize$max_peakwidth + additional
          } else {
            params$no_optimization$max_peakwidth <- params$no_optimization$max_peakwidth + additional
		  }
		}
      }
    }
	
	params <- attachparams(params$to_optimize, params$no_optimization)	    
    iterator <- iterator + 1
                 
  }
  return(history)

}
