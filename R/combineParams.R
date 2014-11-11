combineParams <-
function(params_1, params_2) {
  len <- max(unlist(sapply(params_1, length)))
  #num_params <- length(params_1)
  
  p_names <- c(names(params_1), names(params_2))
  matchedFilter <- "fwhm" %in% p_names
  
  for(i in 1:length(params_2)) {
    new_index <- length(params_1) + 1
	fact <- params_2[[i]]
	params_1[[new_index]] <- fact
	if(matchedFilter) {
	  if(p_names[new_index] == "sigma" && fact == 0) { #update values for sigma if zero
	    if("fwhm" %in% names(params_1)) {
	      params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
	    } else {
	      params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
	    }	
	  } else {
	    if(p_names[new_index] == "mzdiff" && fact == 0) { #update values for mzdiff if zero		
		  if("step" %in% names(params_1)) {
	        if("steps"  %in% names(params_1)) {
	          params_1[[new_index]][1:len] <- 0.8-params_1$step*params_1$steps
	        } else {
	          params_1[[new_index]][1:len] <- 0.8-params_1$step*params_2$steps
	        }	
	      } else {
	        if("steps"  %in% names(params_1)) {
	          params_1[[new_index]][1:len] <- 0.8-params_2$step*params_1$steps
	        } else {
	          params_1[[new_index]][1:len] <- 0.8-params_2$step*params_2$steps
	        }	
	      }		  
		} else {          
          params_1[[new_index]][1:len] <- fact
		}
	  }
	} else {
      params_1[[new_index]][1:len] <- fact
	}
  } 

  names(params_1) <- p_names   
  return(params_1)

}
