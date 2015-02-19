attachList <-
function(params_1, params_2) {
  params <- params_1
  for(factor in params_2)
	params[[length(params)+1]] <- factor
	  
  names(params) <- c(names(params_1), names(params_2))
  return(params)
}


checkParams <-
function(params, quantitative_parameters, qualitative_parameters, unsupported_parameters) { 

  if(length(typeCastParams(params)$to_optimize)==0) {
    stop("No parameters for optimization specified; stopping!")  
  }

  for(i in 1:length(params)) {
	param <- params[[i]]
	name <- names(params)[i]
	if(name %in% unsupported_parameters) {
	  stop(paste("The parameter", name, "is not supported! Please remove from parameters; stopping!"))
	}
	if(name %in% qualitative_parameters) {
	  if(length(param) == 0) {
	    stop(paste("The parameter", name, "has no value set! Please specify; stopping!"))
	  }
	  if(length(param) > 1) {
	    stop(paste("Optimization of parameter", name, "no supported! Please specify only one value; stopping!"))
	  }
	}
    if(name %in% quantitative_parameters) {
	  if(length(param) == 0) {
	    stop(paste("The parameter", name, "has no value set! Please specify between one and two; stopping!"))
	  }
	  if(length(param) > 2) {
	    stop(paste("Too many values for parameter", name, "! Please specify only one or two; stopping!"))
	  }
	}
  } 
  missing_params <- which(!(c(quantitative_parameters, qualitative_parameters) %in% names(params)))
  if(length(missing_params > 0)) {
   stop(paste("The parameter(s)", paste(c(quantitative_parameters, qualitative_parameters)[missing_params], collapse=", "), "are missing! Please specify; stopping!"))
  }
  
}


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


getResponses <-
function(slices, model) {

   if(is.null(slices)) {
     slices= matrix(nrow=1, ncol=0)
   }
   add_col <- array(0, dim=c(nrow(slices), 2))
   slices <- cbind(add_col, slices)
   
   colnames(slices) <- paste("x", 1:ncol(slices), sep="")
   values <- apply(X=slices, MARGIN=1, FUN=getMaxima, model) 
    
   ret <- cbind(t(values), slices[,-c(1,2)]) 
   colnames(ret) <- c("response", paste("x", 1:ncol(slices), sep=""))
    
   return(ret)
}


createModel <-
function(design, params, resp) {
  design$resp <- resp
  formula <- as.formula(paste("resp ~ SO(", paste("x", 1:length(params), sep="", collapse=","), ")", sep="")) 
  return(rsm(formula, data=design)) 
}


decode <-
function(value, bounds) {
  if(is.na(value))
    value <- 1
  x <- (value+1)/2
  x <- (x*(max(bounds)-min(bounds))) + min(bounds)
  
  return(x)
}


decodeAll <-
function(values, params) {

  ret <- rep(0, length(params))
  for(i in 1:length(params))
    ret[i] <- decode(values[i], params[[i]])
  
  names(ret) <- names(params)
  
  return(ret)
}

encode <-
function(value, bounds) {
  x <- (value - min(bounds)) / (max(bounds) - min(bounds))
  
  return(x*2-1)
}


getBbdParameter <-
function(params) {
 
  lower_bounds <- unlist(lapply(X=params, FUN=function(x) x[1]))
  higher_bounds <- unlist(lapply(X=params, FUN=function(x) x[2]))
  
  steps <- (higher_bounds - lower_bounds)/2

  x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", (lower_bounds + steps), ")/", steps, sep="")
  formulae <- list()
  for(i in 1:length(x))
    formulae[[i]] <- as.formula(x[i])  
  
  design <- bbd(length(params), n0 = 1, randomize = FALSE, coding = formulae)

  return(design)
  
}


getCcdParameter <-
function(params) {
 
  lower_bounds <- unlist(lapply(X=params, FUN=function(x) x[1]))
  higher_bounds <- unlist(lapply(X=params, FUN=function(x) x[2]))
  
  steps <- (higher_bounds - lower_bounds)/2
  
  x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", (lower_bounds + steps), ")/", steps, sep="")
  formulae <- list()
  for(i in 1:length(x))
    formulae[[i]] <- as.formula(x[i])  
  
  design <- ccd(length(params), n0 = 1, alpha = "face", randomize = FALSE, inscribed = TRUE, coding = formulae)

  return(design)
  
}


getMaxima <-
function(slice, model) {
  values <- contour(model, ~ x1 * x2, at=slice, plot.it=FALSE, decode=FALSE) 
    
  z_values <- values[[1]]$z
  max_z <- max(z_values)
  maxima <- which(z_values==max_z, arr.ind=T)
  ret <- c(max_z, 0,0)
  if(length(maxima)==2){
    ret[2] <- values[[1]]$x[maxima[1]]
    ret[3] <- values[[1]]$y[maxima[2]]     
  } else {
    for(i in 1:2) {
      col_max <- maxima[,i]
	    if(sum(col_max==1) > 0 && sum(col_max==length(values[[1]][[i]])) > 0) 
	      ret[(i+1)] <- NA
      else
   	    ret[(i+1)] <- values[[1]][[i]][max(col_max)]
	  }
  }
  
  return(ret)
}


getMaximumExperiment <-
function(model) {  
  dimensions <- length(model$coding)   
  slices <- getSlices(dimensions-2)
  mat <- getResponses(slices, model)  
  return(mat[which.max(mat[,1]),])
}


getSlices <-
function(dimensions, slice=NULL) {
  ret <- c()
  if(length(slice) == dimensions) {
    return(slice)
  } else {
    values <- seq(-1,1,0.2)
    for(val in values)
      ret <- rbind(ret, getSlices(dimensions, c(slice, val)))
  }
  
  return(ret)
}


plotContours <-
function(model, maximum_slice, plot_name) {

  plots <- c()
  for(i in 1:(length(maximum_slice)-1)) {
    for(j in (i+1):length(maximum_slice)) {
      plots <- c(plots, as.formula(paste("~ x", i, "* x", j, sep="")))
    } 
  }
    
  plot_rows <- round(sqrt(length(plots)))
  plot_cols <- if(plot_rows==1){length(plots)}else{ceiling(sqrt(length(plots)))}

  if(!is.null(plot_name)) {
    plot_name = paste(plot_name, ".jpg", sep="")
    jpeg(plot_name, width=4*plot_cols, height=2*plot_rows+2, units="in", res=c(200,200))
  } else 
    dev.new(width=4*plot_cols, height=2*plot_rows+2)
  par(mfrow=c(plot_rows, plot_cols), oma=c(3,0,2,0))  
  contours <- contour(model, plots, image=TRUE, at=maximum_slice)
  
  if(!is.null(plot_name)) {
    dev.off()
  }

}


toMatrix <-
function(data) {

  if(!is.matrix(data)) {
    tmp <- names(data)
    data <- matrix(data, nrow=1)
    colnames(data) <- tmp  
  } 
  
  return(data)  
  
}


## Removed after switching away from MPI
## startSlaves <-
## function(nSlaves) {
##   mpi.spawn.Rslaves(nslaves=nSlaves)
                                                                                
##   .Last <- function() {
##     if (is.loaded("mpi_initialize")) {
##       if (mpi.comm.size(1) > 0) {
## 	    print("Please use mpi.close.Rslaves() to close slaves.")
## 	    mpi.close.Rslaves()
## 	  }
##       print("Please use mpi.quit() to quit R")
##       .Call("mpi_finalize", PACKAGE="Rmpi")
##     }
##   }
## }


typeCastParams <-
function(params) {
  ret_1 <- list()
  ret_2 <- list()  
  ret <- list()
  for(i in  1:length(params)) {
    factor <- params[[i]]
    if(length(factor) == 2) {
	  ret_1[[(length(ret_1)+1)]] <- factor
	  names(ret_1)[length(ret_1)] <- names(params)[i]
	} else {	
	  ret_2[[(length(ret_2)+1)]] <- factor
	  names(ret_2)[length(ret_2)] <- names(params)[i]
	}
  }	
  ret$to_optimize <- ret_1
  ret$no_optimization <- ret_2
  
  return(ret)
}


writeRScript <-
function(peakPickingSettings, retCorGroupSettings, nSlaves) {
  cat("library(xcms)\n")
  cat("library(Rmpi)\n")

  if(is.null(peakPickingSettings$step)) {     #centWave     		
    cat(paste("xset <- xcmsSet(method=\"centWave\", peakwidth=c(", 
              peakPickingSettings$min_peakwidth, ", ", peakPickingSettings$max_peakwidth,
              "), ppm=", peakPickingSettings$ppm, ", noise=", peakPickingSettings$noise, 
			  ", snthresh=", peakPickingSettings$snthresh, ", mzdiff=", peakPickingSettings$mzdiff,
			  ", prefilter=c(", peakPickingSettings$prefilter, ", ", peakPickingSettings$value_of_prefilter,
			  "), mzCenterFun=\"", peakPickingSettings$mzCenterFun, "\", integrate=", peakPickingSettings$integrate,
			  ", fitgauss=", peakPickingSettings$fitgauss, ", verbose.columns=", peakPickingSettings$verbose.columns,
			  ", nSlaves=", nSlaves, ")\n", sep=""))
                  
  } else { #matchedFilter  
    cat(paste("xset <- xcmsSet(method=\"matchedFilter\", fwhm=", 
	         peakPickingSettings$fwhm, ", snthresh=",peakPickingSettings$snthresh,
             ", step=", peakPickingSettings$step, ", steps=", peakPickingSettings$steps,
             ", sigma=", peakPickingSettings$sigma, ", max=", peakPickingSettings$max, 
             ", mzdiff=", peakPickingSettings$mzdiff, ", index=", peakPickingSettings$index,
	         ")\n", sep=""))   
  }
	  
	  
  cat(paste("xset <- retcor(xset, method=\"obiwarp\", plottype=\"", retCorGroupSettings$plottype, 
            "\", distFunc=\"", retCorGroupSettings$distFunc, "\", profStep=", retCorGroupSettings$profStep, 
			", center=", retCorGroupSettings$center, ", response=", retCorGroupSettings$response, 
		    ", gapInit=", retCorGroupSettings$gapInit, ", gapExtend=", retCorGroupSettings$gapExtend,
		    ", factorDiag=", retCorGroupSettings$factorDiag, ", factorGap=", retCorGroupSettings$factorGap, 
			", localAlignment=", retCorGroupSettings$localAlignment, ")\n", sep=""))
  	   
  cat(paste("xset <- group(xset, method=\"density\", bw=", retCorGroupSettings$bw, 
            ", mzwid=", retCorGroupSettings$mzwid, ", minfrac=", retCorGroupSettings$minfrac, 
            ", minsamp=", retCorGroupSettings$minsamp, ", max=", retCorGroupSettings$max,
	        ")\n", sep=""))	 
	  
  cat(paste("xset <- fillPeaks(xset, nSlaves=", nSlaves, ")\n", sep=""))	
	
}


writeParamsTable <- 
function(peakPickingSettings, retCorGroupSettings, file, ...) {
  write.table(combineParams(peakPickingSettings, retCorGroupSettings), file, ...)
}