attachparams <-
function(params_1, params_2) {
  params <- params_1
  for(factor in params_2)
	params[[length(params)+1]] <- factor
	  
  names(params) <- c(names(params_1), names(params_2))
  return(params)
}
