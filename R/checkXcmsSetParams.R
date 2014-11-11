checkXcmsSetParams <-
function(params) {
  if(is.null(params$step)) {     #centWave   
    quantitative_parameters <- c("ppm", "min_peakwidth", "max_peakwidth", "snthresh", "mzdiff", "noise", "prefilter", "value_of_prefilter")
    qualitative_parameters <- c("integrate", "fitgauss", "verbose.columns", "mzCenterFun")
    unsupported_parameters <- c("scanrange", "sleep", "ROI.list")
  } else {    
    quantitative_parameters <- c("fwhm", "sigma", "max", "snthresh", "step", "steps", "mzdiff")
    qualitative_parameters <- c("index")
    unsupported_parameters <- c("sleep")  
  } 
  checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
}
