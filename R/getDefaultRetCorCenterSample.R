getDefaultRetCorCenterSample <-
function(xset) {
  ret <- NULL
  for(i in 1:length(xset@filepaths)) {
    ret <- c(ret, sum(xset@peaks[,"sample"] == i))
  }
  return(which.max(ret))
}
