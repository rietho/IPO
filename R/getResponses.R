getResponses <-
function(slices, model) {
   add_col <- array(0, dim=c(nrow(slices), 2))
   slices <- cbind(add_col, slices)
   
   colnames(slices) <- paste("x", 1:ncol(slices), sep="")
   values <- apply(X=slices, MARGIN=1, FUN=getMaxima, model) 
    
   ret <- cbind(t(values), slices[,-c(1,2)]) 
   colnames(ret) <- c("response", paste("x", 1:ncol(slices), sep=""))
    
   return(ret)
}
