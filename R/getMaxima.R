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
