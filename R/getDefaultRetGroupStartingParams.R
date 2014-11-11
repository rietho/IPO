getDefaultRetGroupStartingParams <-
function(distfunc="cor_opt", high_resolution=TRUE) {

  ret <- NULL
  if(!is.null(distfunc)) {
    if(distfunc=="cor")
      ret <- (list(distFunc="cor", gapInit=c(0.0, 0.4), gapExtend=c(2.1, 2.7)))
	if(distfunc=="cor_opt")
	  ret <- (list(distFunc="cor_opt", gapInit=c(0.0, 0.4), gapExtend=c(2.1, 2.7)))
    if(distfunc=="cov")
	  ret <- (list(distFunc="cov", gapInit=c(0.0, 0.4), gapExtend=c(11.4, 12.0)))
	if(distfunc=="prd")
	  ret <- (list(distFunc="prd", gapInit=c(0.0, 0.4), gapExtend=c(7.5, 8.1)))
	if(distfunc=="euc")
	  ret <- (list(distFunc="euc", gapInit=c(0.7, 1.1), gapExtend=c(1.5, 2.1)))

    ret$profStep <- c(0.7, 1)
    ret$plottype <- "none"
    ret$response <- 1
    ret$factorDiag <- 2
    ret$factorGap <- 1
    ret$localAlignment <- 0
    #ret$initPenalty <- 0
	  
  } else {
	ret <- list()
  }

	#grouping parameter
  ret$bw <- c(22,38)
  ret$minfrac <- c(0.3, 0.7)
  ret$mzwid <- c(0.15, 0.35)
  if(high_resolution)
    ret$mzwid <- c(0.015, 0.035)
  ret$minsamp <- 1
  ret$max <- 50
  #ret$sleep <- 0
	  
  return(ret)

}
