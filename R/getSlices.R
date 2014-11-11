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
