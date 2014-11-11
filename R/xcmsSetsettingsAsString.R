xcmsSetsettingsAsString <-
function(parameters) {  
  
  return(sprintf("xcmsSet(method='centWave', ppm=%s, peakwidth=c(%s,%s), snthresh=%s, prefilter=c(%s,%s), mzCenterFun='%s', integrate=%s, mzdiff=%s, fitgauss=%s, noise=%s, verbose.columns=%s)", 
          prettyNum(parameters$ppm), prettyNum(parameters$min_peakwidth), 
          prettyNum(parameters$max_peakwidth), prettyNum(parameters$snthresh), 
          prettyNum(parameters$prefilter), prettyNum(parameters$value_of_prefilter), 
          parameters$mzCenterFun, prettyNum(parameters$integrate), 
          prettyNum(parameters$mzdiff), parameters$fitgauss, 
          prettyNum(parameters$noise), parameters$verbose.columns))
}
