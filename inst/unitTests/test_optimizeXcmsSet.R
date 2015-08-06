test_ipo <- function() {
  mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
  
  #checking peak picking optimization 
  paramsPP <- getDefaultXcmsSetStartingParams()
  paramsPP$mzdiff <- -0.001
  paramsPP$min_peakwidth <- c(7,14)
  paramsPP$max_peakwidth <- c(20,30)
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir="microtofq_IPO", 
                              checkPeakShape="borderIntensity")
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 222, 69, 77)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 85.92754, tolerance=1e-3) 
  
  
  #checking findIsotopes.IPO
  xset <- resultPP$best_settings$xset
  PPS <- calcPPS(xset, "IPO", checkPeakShape="none")
  checkTrue(all(PPS[1:4]== c(0, 222, 67, 81)))
  checkEqualsNumeric(PPS[5], 97.92537, tolerance=1e-3) 
  PPS <- calcPPS(xset, "IPO", checkPeakShape="sinusCurve")
  checkTrue(all(PPS[1:4]== c(0, 222, 123, 10)))
  checkEqualsNumeric(PPS[5], 0.8130081, tolerance=1e-3) 
  PPS <- calcPPS(xset, "IPO", checkPeakShape="normalDistr")
  checkTrue(all(PPS[1:4]== c(0, 222, 125, 8)))
  checkEqualsNumeric(PPS[5], 0.512, tolerance=1e-3) 
  
  
  #checking peak picking optimization using CAMERA isotope identification
  resultPPCamera <- optimizeXcmsSet(mzmlfile, paramsPP, isotopeIdentification="CAMERA", 
                                    subdir="microtofq_CAMERA", ppm=15, maxcharge=2)
  checkTrue(all(resultPPCamera$best_settings$result[1:4]== c(0, 162, 52, 75)))
  checkEqualsNumeric(resultPPCamera$best_settings$result[5], 108.1731, tolerance=1e-3)
  
  
  #checking single parameter peak picking optimization
  paramsPP$max_peakwidth <- 15
  paramsPP$ppm <- 50
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir=NULL, 
                              checkPeakShape="borderIntensity") 
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 218, 83, 68)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 55.71084, tolerance=1e-3)  
  
  #checking retention time correction and grouping optimization
  mtbls2files <- list.files(paste(find.package("mtbls2"), "/mzData", sep=""), 
                            full.names=TRUE)
  xset <- xcmsSet(mtbls2files[1:3], method="centWave", peakwidth=(c(12, 30)), ppm=30)
  
  
  #checking obiwarp
  paramsRG <- getDefaultRetGroupStartingParams()
  paramsRG$minfrac <- 1
  paramsRG$profStep <- 1
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir="mtbls2Obiwarp")

  TV <- resultRG[[4]]$target_value
  checkTrue(all(unlist(TV)[-c(5)]== c(1, 2656, 1, 7054336, 1)))
  checkEqualsNumeric(unlist(TV)[5], 1593.865, tolerance=1e-2)

  
  #checking single parameter retention time correction and grouping optimization
  paramsRG$gapInit <- 0.4
  paramsRG$gapExtend <- 2.4
  paramsRG$mzwid <- 0.01 
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir="mtbls2")
  
  TV <- resultRG[[4]]$target_value
  checkTrue(all(unlist(TV)[-c(5)]== c(1, 2652, 1, 7033104, 1)))
  checkEqualsNumeric(unlist(TV)[5], 1586.247, tolerance=1e-2)
  
  #checking loess
  paramsRG <- getDefaultRetGroupStartingParams("loess")
  paramsRG$extra <- 0
  paramsRG$missing <- 0
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir="mtbls2Loess")
  
  TV <- resultRG[[3]]$target_value
  checkTrue(all(unlist(TV)[-c(5)]== c(1, 2664, 1, 7096896, 1)))
  checkEqualsNumeric(unlist(TV)[5], 1547.593, tolerance=3e-2)
  
}
