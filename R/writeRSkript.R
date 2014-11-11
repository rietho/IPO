writeRSkript <-
function(peakPickingSettings, retCorGroupSettings, nSlaves) {
  cat("library(xcms)\n")
  cat("library(Rmpi)\n")

  if(is.null(peakPickingSettings$step)) {     #centWave     		
    cat(paste("xset <- xcmsSet(method=\"centWave\", peakwidth=c(", 
              peakPickingSettings$min_peakwidth, ", ", peakPickingSettings$max_peakwidth,
              "), ppm=", peakPickingSettings$ppm, ", noise=", peakPickingSettings$noise, 
			  ", snthresh=", peakPickingSettings$snthresh, ", mzdiff=", peakPickingSettings$mzdiff,
			  ", prefilter=c(", peakPickingSettings$prefilter, ", ", peakPickingSettings$value_of_prefilter,
			  "), mzCenterFun=\"", peakPickingSettings$mzCenterFun, "\", integrate=", peakPickingSettings$integrate,
			  ", fitgauss=", peakPickingSettings$fitgauss, ", verbose.columns=", peakPickingSettings$verbose.columns,
			  ", nSlaves=", nSlaves, ")\n", sep=""))
                  
  } else { #matchedFilter  
    cat(paste("xset <- xcmsSet(method=\"matchedFilter\", fwhm=", 
	         peakPickingSettings$fwhm, ", snthresh=",peakPickingSettings$snthresh,
             ", step=", peakPickingSettings$step, ", steps=", peakPickingSettings$steps,
             ", sigma=", peakPickingSettings$sigma, ", max=", peakPickingSettings$max, 
             ", mzdiff=", peakPickingSettings$mzdiff, ", index=", peakPickingSettings$index,
	         ")\n", sep=""))   
  }
	  
	  
  cat(paste("xset <- retcor(xset, method=\"obiwarp\", plottype=\"", retCorGroupSettings$plottype, 
            "\", distFunc=\"", retCorGroupSettings$distFunc, "\", profStep=", retCorGroupSettings$profStep, 
			", center=", retCorGroupSettings$center, ", response=", retCorGroupSettings$response, 
		    ", gapInit=", retCorGroupSettings$gapInit, ", gapExtend=", retCorGroupSettings$gapExtend,
		    ", factorDiag=", retCorGroupSettings$factorDiag, ", factorGap=", retCorGroupSettings$factorGap, 
			", localAlignment=", retCorGroupSettings$localAlignment, ")\n", sep=""))
  	   
  cat(paste("xset <- group(xset, method=\"density\", bw=", retCorGroupSettings$bw, 
            ", mzwid=", retCorGroupSettings$mzwid, ", minfrac=", retCorGroupSettings$minfrac, 
            ", minsamp=", retCorGroupSettings$minsamp, ", max=", retCorGroupSettings$max,
	        ")\n", sep=""))	 
	  
  cat(paste("xset <- fillPeaks(xset, nSlaves=", nSlaves, ")\n", sep=""))	
	
}
