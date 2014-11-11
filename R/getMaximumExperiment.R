getMaximumExperiment <-
function(model) {  
  dimensions <- length(model$coding)   
  slices <- getSlices(dimensions-2)
  mat <- getResponses(slices, model)  
  return(mat[which.max(mat[,1]),])
}
