typeCastFactor <-
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
