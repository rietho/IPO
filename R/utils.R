attachList <-
function(params_1, params_2) {
  params <- params_1
  for(factor in params_2)
	params[[length(params)+1]] <- factor
	  
  names(params) <- c(names(params_1), names(params_2))
  return(params)
}


checkParams <-
function(params, quantitative_parameters, qualitative_parameters, 
         unsupported_parameters) { 

  if(length(typeCastParams(params)$to_optimize)==0) {
    stop("No parameters for optimization specified; stopping!")  
  }

  for(i in 1:length(params)) {
	  param <- params[[i]]
	  name <- names(params)[i]
	  if(name %in% unsupported_parameters) {
	    stop(paste("The parameter", name, "is not supported! Please remove
	               from parameters; stopping!"))
	  }
	  if(name %in% qualitative_parameters) {
	    if(length(param) == 0) {
	      stop(paste("The parameter", name, "has no value set!
	                 Please specify; stopping!"))
	    }
	    if(length(param) > 1) {
	      stop(paste("Optimization of parameter", name, "not supported!
	                 Please specify only one value; stopping!"))
	    }
	  }
      if(name %in% quantitative_parameters) {
	      if(length(param) == 0) {
	        stop(paste("The parameter", name, "has no value set!
	                   Please specify between one and two; stopping!"))
	      } 
	      if(length(param) > 2) {
	        stop(paste("Too many values for parameter", name, "!
	                   Please specify only one or two; stopping!"))
	      }
	    }
  } 
  missing_params <- which(!(c(quantitative_parameters, 
                              qualitative_parameters) %in% names(params)))
  if(length(missing_params > 0)) {
    stop(paste("The parameter(s)", paste(c(quantitative_parameters, qualitative_parameters)
                                         [missing_params], collapse=", "), 
               "is/are missing! Please specify; stopping!"))
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
	    #update values for sigma if zero
	    if(p_names[new_index] == "sigma" && fact == 0) {
	      if("fwhm" %in% names(params_1)) {
	        params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
	      } else {
	        params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
	      }	
	    } else {
	      #update values for mzdiff if zero	
	      if(p_names[new_index] == "mzdiff" && fact == 0) { 	
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

createModel <-
function(design, params, resp) {
  design$resp <- resp
  if(length(params) > 1) {
    formula <- as.formula(paste("resp ~ SO(", paste("x", 1:length(params), sep="", 
                                                    collapse=","), ")", sep="")) 
    model <- rsm(formula, data=design) 
  } else {
    param_name <- names(params)[1]
    formula <- as.formula(paste("resp ~ ", param_name, " + ", 
                                param_name, " ^ 2", sep="")) 
    model <- lm(formula, data=design) 
    model$coding <- list(x1=as.formula(paste(param_name, "~ x1"))) 
    names(model$coding) <- param_name
    #attr(model, "class") <- c("rsm", "lm")
  }
  return(model)  
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

  x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", 
             (lower_bounds + steps), ")/", steps, sep="")
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
  
  x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", 
             (lower_bounds + steps), ")/", steps, sep="")
  formulae <- list()
  for(i in 1:length(x))
    formulae[[i]] <- as.formula(x[i])  
  
  design <- ccd(length(params), n0 = 1, alpha = "face", randomize = FALSE, 
                inscribed = TRUE, coding = formulae)

  return(design)
  
}



getMaximumLevels <-
function(model) {  
  dimensions <- length(model$coding)   
  #slices <- getSlices(dimensions-2)
  #mat <- getResponses(slices, model)
  #testdata <- getTestData(dimensions)
  #if(dimensions==1)
  #  names(testdata) <- names(model$coding)
  #
  #return(getMaxSettings(testdata, model))
  
  if(dimensions > 6) {
    testSpace <- seq(-1,1,0.2)
  } else { 
    testSpace <- seq(-1,1,0.1)
  }
  
  testSize <- 1000000
  testAmount <- length(testSpace)^dimensions
  i <- 1
  max <- rep(-1, dimensions+1)
  while(i < testAmount) {
    testdata <- expand.grid.subset(i:(i+testSize), testSpace, dimensions)
    if(dimensions==1)
      names(testdata) <- names(model$coding)
    max_tmp <- getMaxSettings(testdata, model)
    if(max_tmp[1]>max[1])
      max <- max_tmp
    i <- i + testSize + 1
  }
  
  return(max)
  
}

getMaxSettings <-  function(testdata, model) {
    
  response <- predict(model, testdata)
  max_response <- max(response)
  max_settings <- testdata[response==max_response,,drop=FALSE]
  ret <- max_response
  for(i in 1:ncol(testdata)) {
    levels <- max_settings[,i]
    if(sum(c(-1,1) %in% levels)==2)
       ret <- cbind(ret, NA)
    else
      ret <- cbind(ret,levels[1])
  }
    
  colnames(ret) <- c("response", paste("x", 1:ncol(testdata), sep=""))
    
  return(ret)
}


expand.grid.subset  <- function(subset, sequence, dimensions) { 
  
  vars <- list()
  for(i in 1:dimensions) {
    vars[[i]] <- sequence
  }
  names(vars) <- paste("x", 1:dimensions, sep="")
  
  maximumSubset <- length(sequence)^dimensions
  subset <- min(maximumSubset,min(subset)):min(maximumSubset, max(subset))
  
  #nv <-  #length(vars) 
  lims <- sapply(vars,length) 
  stopifnot(length(lims) > 0, subset <= prod(lims), length(names(vars)) == dimensions) 
  res <- structure(vector("list",dimensions), .Names = names(vars)) 
  if (dimensions > 1) {
    for(i in dimensions:2) { 
      f <- prod(lims[1:(i-1)]) 
      res[[i]] <- vars[[i]][(subset - 1)%/%f + 1] 
      subset <- (subset - 1)%%f + 1 
    } 
  }
  res[[1]] <- vars[[1]][subset] 
  as.data.frame(res) 
} 


#getTestData <- function(parameters) {
#  step=0.1
#  if(parameters < 2)
#    step <- 0.05
#  if(parameters > 5)
#    step <- 0.2

#  m <- matrix(rep(seq(-1,1,step),parameters), byrow=FALSE, ncol=parameters,   dimnames=list(NULL, paste("x", 1:parameters,  sep=""))) 
  #colnames(m) <-  paste("x", 1:parameters,  sep="")
  
#  return(expand.grid(data.frame(m)))
#}


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
    jpeg(plot_name, width=4*plot_cols, height=2*plot_rows+2, 
         units="in", res=c(200,200))
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
  message("library(xcms)\n")
  message("library(Rmpi)\n")

  if(is.null(peakPickingSettings$step)) {     #centWave     		
    message(paste("xset <- xcmsSet(method=\"centWave\", 
                  peakwidth=c(", peakPickingSettings$min_peakwidth, ", ", 
                  peakPickingSettings$max_peakwidth,
                  "), ppm=", peakPickingSettings$ppm, 
                  ", noise=", peakPickingSettings$noise, 
                  ", snthresh=", peakPickingSettings$snthresh, 
                  ", mzdiff=", peakPickingSettings$mzdiff,
                  ", prefilter=c(", peakPickingSettings$prefilter, 
                  ", ", peakPickingSettings$value_of_prefilter,			  
                  "), mzCenterFun=\"", peakPickingSettings$mzCenterFun, 
                  "\", integrate=", peakPickingSettings$integrate,
                  ", fitgauss=", peakPickingSettings$fitgauss,
                  ", verbose.columns=", peakPickingSettings$verbose.columns,
                  ", nSlaves=", nSlaves, ")", sep=""))
                  
  } else { #matchedFilter  
    message(paste("xset <- xcmsSet(method=\"matchedFilter\", fwhm=", 
                  peakPickingSettings$fwhm, 
                  ", snthresh=",peakPickingSettings$snthresh,
                  ", step=", peakPickingSettings$step, 
                  ", steps=", round(peakPickingSettings$steps),
                  ", sigma=", peakPickingSettings$sigma, 
                  ", max=", round(peakPickingSettings$max), 
                  ", mzdiff=", peakPickingSettings$mzdiff,
                  ", index=", peakPickingSettings$index,
                  ", nSlaves=", nSlaves, ")", sep=""))   
  }
	  
  if(retCorGroupSettings$retcorMethod == "loess")	{
    
    message(paste("xset <- group(xset, method=\"density\", bw=", 
                  retCorGroupSettings$bw, 
                  ", mzwid=", retCorGroupSettings$mzwid, 
                  ", minfrac=", retCorGroupSettings$minfrac,
                  ", minsamp=", round(retCorGroupSettings$minsamp), 
                  ", max=", round(retCorGroupSettings$max), ")", sep=""))	 
    
    message(paste("xset <- retcor(xset", 
                  ", missing=", round(retCorGroupSettings$missing), 
                  ", extra=", round(retCorGroupSettings$extra),
                  ", span=", retCorGroupSettings$span, 
                  ", smooth=\"", retCorGroupSettings$smooth, 
                  "\", family=\"", retCorGroupSettings$family,         
                  "\", plottype=\"", retCorGroupSettings$plottype,
                  "\")", sep=""))	 
  }  
  
  if(retCorGroupSettings$retcorMethod == "obiwarp") {
    message(paste("xset <- retcor(xset, method=\"obiwarp\",
                  plottype=\"", retCorGroupSettings$plottype, 
                  "\", distFunc=\"", retCorGroupSettings$distFunc, 
                  "\", profStep=", retCorGroupSettings$profStep, 
                  ", center=", retCorGroupSettings$center, 
                  ", response=", retCorGroupSettings$response,
                  ", gapInit=", retCorGroupSettings$gapInit,
                  ", gapExtend=", retCorGroupSettings$gapExtend,
                  ", factorDiag=", retCorGroupSettings$factorDiag, 
                  ", factorGap=", retCorGroupSettings$factorGap, 
                  ", localAlignment=", retCorGroupSettings$localAlignment, 
                  ")", sep=""))
  }
  	   
  message(paste("xset <- group(xset, method=\"density\", 
                bw=", retCorGroupSettings$bw, 
                ", mzwid=", retCorGroupSettings$mzwid,
                ", minfrac=", retCorGroupSettings$minfrac, 
                ", minsamp=", retCorGroupSettings$minsamp,
                ", max=", retCorGroupSettings$max, ")\n", sep=""))	 
	  
  message(paste("xset <- fillPeaks(xset, nSlaves=", nSlaves, ")", sep=""))	
	
}


writeParamsTable <- 
function(peakPickingSettings, retCorGroupSettings, file, ...) {
  write.table(combineParams(peakPickingSettings, retCorGroupSettings), file, ...)
}