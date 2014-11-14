optimizeRetGroup <-
function(xset, params=getDefaultRetGroupStartingParams(), nSlaves=4, subdir="IPO") {
                                                 
  library(xcms)
  iterator = 1 
  history <- list()  
  best_range <- 0.25

  if(!file.exists(subdir))
    dir.create(file.path(getwd(), subdir))
	
  if(is.null(params$center))
    params$center <- getDefaultRetCorCenterSample(xset) 

  checkRetGroupSetParams(params)	
    
  while(iterator < 50) {
    cat("\n")
    cat("\n")
    cat("\n")
    cat("starting new DoE with:\n")
    print(params)
        
    retcor_result <- retGroupCalcExperiments(params, xset, nSlaves)  
                       
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
