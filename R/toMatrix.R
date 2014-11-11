toMatrix <-
function(data) {

  if(!is.matrix(data)) {
    tmp <- names(data)
    data <- matrix(data, nrow=1)
    colnames(data) <- tmp  
  } 
  
  return(data)  
  
}
