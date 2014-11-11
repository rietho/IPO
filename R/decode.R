decode <-
function(value, bounds) {
  if(is.na(value))
    value <- 1
  x <- (value+1)/2
  x <- (x*(max(bounds)-min(bounds))) + min(bounds)
  
  return(x)
}
