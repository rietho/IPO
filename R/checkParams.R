checkParams <-
function(params, quantitative_parameters, qualitative_parameters, unsupported_parameters) { 
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
