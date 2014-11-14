getRGTVValues <-
function(xset, exp_index=1, retcor_penalty=1) {
  features <- nrow(xset@groups)

  good_groups <- sum(unlist(lapply(X=xset@groupidx, FUN=function(x, xset) {ifelse(length(unique(xset@peaks[x,"sample"]))==length(xset@filepaths) & 
                                   length(xset@peaks[x,"sample"])==length(xset@filepaths),1,0)}, xset)))
  bad_groups <- nrow(xset@groups) - good_groups

  relative_rt_diff <- c()
  
  if(nrow(xset@groups) > 0) {
    for(i in 1:nrow(xset@groups)) {
      feature_rtmed <- xset@groups[i, "rtmed"]
	    relative_rt_diff <- c(relative_rt_diff, mean(abs(feature_rtmed - xset@peaks[xset@groupidx[[i]], "rt"])/feature_rtmed))
    }
  } else {
    relative_rt_diff <- 1
  }
  
  ARTS <- (mean(relative_rt_diff)) * retcor_penalty
  
  ret <- list(exp_index=exp_index, good_groups=good_groups, bad_groups=bad_groups, mean_rel_rt_diff=ARTS)
  
  ret$retcor_done = retcor_penalty        
  
  return(ret)  
}
