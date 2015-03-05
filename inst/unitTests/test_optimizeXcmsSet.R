test_ipo <- function() {
  mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
  
  paramsPP <- getDefaultXcmsSetStartingParams()
  paramsPP$mzdiff <- -0.001
  paramsPP$min_peakwidth <- c(7,14)
  paramsPP$max_peakwidth <- c(20,30)
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir="microtofq_IPO", checkBorderIntensity=TRUE)
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 225, 64, 82)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 105.0625, tolerance=0.000001)
  
  
  resultPPCamera <- optimizeXcmsSet(mzmlfile, paramsPP, isotopeIdentification="CAMERA", subdir="microtofq_CAMERA", ppm=15, maxcharge=2)
  
  checkTrue(all(resultPPCamera$best_settings$result[1:4]== c(0, 216, 64, 96)))
  checkEqualsNumeric(resultPPCamera$best_settings$result[5], 144, tolerance=0.000001)
  
}
