getDefaultXcmsSetStartingParams <-
function(method="centWave") {
  if(method=="centWave")
    return(list(min_peakwidth=c(10,30), max_peakwidth=c(35,65), ppm=c(15,35),
              mzdiff=c(-0.001, 0.01), snthresh=10, noise=0, prefilter=3, 
			        value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
			        fitgauss=FALSE, verbose.columns=FALSE))
			  
  if(method=="matchedFilter")
    return(list(fwhm=c(25,35), snthresh=c(3,17), step=c(0.05, 0.15), steps=c(1,3), 
           sigma=0, max=5, mzdiff=0, index=FALSE)) 
  
}
