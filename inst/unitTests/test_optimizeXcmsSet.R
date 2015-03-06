test_ipo <- function() {
  mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
  
  #checking peak picking optimization 
  paramsPP <- getDefaultXcmsSetStartingParams()
  paramsPP$mzdiff <- -0.001
  paramsPP$min_peakwidth <- c(7,14)
  paramsPP$max_peakwidth <- c(20,30)
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir="microtofq_IPO", checkBorderIntensity=TRUE)
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 222, 69, 77)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 85.92754, tolerance=0.000001)
    
  
  #checking peak picking optimization using CAMERA isotope identification
  resultPPCamera <- optimizeXcmsSet(mzmlfile, paramsPP, isotopeIdentification="CAMERA", subdir="microtofq_CAMERA", ppm=15, maxcharge=2)
  checkTrue(all(resultPPCamera$best_settings$result[1:4]== c(0, 182, 61, 75)))
  checkEqualsNumeric(resultPPCamera$best_settings$result[5], 92.21311, tolerance=0.000001)
  
  
  #checking single parameter peak picking optimization
  paramsPP$max_peakwidth <- 15
  paramsPP$ppm <- 50
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir=NULL, checkBorderIntensity=FALSE) 
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 218, 78, 78)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 78, tolerance=0.000001)
  
  
  
  #checking retention time correction and grouping optimization
  mtbls2files <- list.files(paste(find.package("mtbls2"), "/mzdata", sep=""), full.names=TRUE)
  xset <- xcmsSet(mtbls2files[1:3], method="centWave", peakwidth=(c(12, 30)), ppm=30)
  
  paramsRG <- getDefaultRetGroupStartingParams()
  paramsRG$minfrac <- 1
  paramsRG$profStep <- 1
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir="mtbls2")

  TV <- resultRG[[4]]$target_value
  checkTrue(all(unlist(TV)[-4]== c(1, 2656, 1, 1)))
  checkEqualsNumeric(unlist(TV)[4], 0.0006312929, tolerance=0.000001)
  
  #checking single parameter retention time correction and grouping optimization
  paramsRG$gapInit <- 0.4
  paramsRG$gapExtend <- 2.4
  paramsRG$mzwid <- 0.01 
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir="mtbls2")
  
  TV <- resultRG[[4]]$target_value
  checkTrue(all(unlist(TV)[-4]== c(1, 2652, 1, 1)))
  checkEqualsNumeric(unlist(TV)[4], 0.0006276349, tolerance=0.000001)
  
}
