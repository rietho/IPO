test_ipo <- function() {
  mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
  
  paramsPP <- getDefaultXcmsSetStartingParams()
  paramsPP$mzdiff <- -0.001
  paramsPP$min_peakwidth <- c(7,14)
  paramsPP$max_peakwidth <- c(20,30)
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir="microtofq_IPO")
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 221, 65, 81)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 100.9385, tolerance=0.000001)
  
  
  resultPPCamera <- optimizeXcmsSet(mzmlfile, paramsPP, isotopeIdentification="CAMERA", subdir="microtofq_CAMERA", ppm=15, maxcharge=2)
  
  checkTrue(all(resultPPCamera$best_settings$result[1:4]== c(0, 215, 62, 97)))
  checkEqualsNumeric(resultPPCamera$best_settings$result[5], 151.7581, tolerance=0.000001)
  
}
