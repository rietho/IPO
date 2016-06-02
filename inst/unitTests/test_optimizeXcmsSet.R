test_ipo <- function() {
  mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
  
  #checking peak picking optimization 
  paramsPP <- getDefaultXcmsSetStartingParams()
  paramsPP$mzdiff <- -0.001
  paramsPP$ppm <- 56
  paramsPP$min_peakwidth <- c(3, 9.5)
  paramsPP$max_peakwidth <- c(10,20)
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir=NULL, nSlaves=2, 
                              checkPeakShape="borderIntensity")
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 225, 64, 82)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 105.0625, tolerance=1e-3) 
  
  
  #checking findIsotopes.IPO
  xset <- resultPP$best_settings$xset
  PPS <- calcPPS(xset, "IPO", checkPeakShape="none")
  checkTrue(all(PPS[1:4]== c(0, 225, 64, 82)))
  checkEqualsNumeric(PPS[5], 105.0625, tolerance=1e-3) 
  PPS <- calcPPS(xset, "IPO", checkPeakShape="sinusCurve")
  checkTrue(all(PPS[1:4]== c(0, 225, 121, 14)))
  checkEqualsNumeric(PPS[5], 1.619835, tolerance=1e-3) 
  PPS <- calcPPS(xset, "IPO", checkPeakShape="normalDistr")
  checkTrue(all(PPS[1:4]== c(0, 225, 123, 12)))
  checkEqualsNumeric(PPS[5], 1.170732, tolerance=1e-3) 
  
  
  #checking peak picking optimization using CAMERA isotope identification
  resultPPCamera <- optimizeXcmsSet(mzmlfile, paramsPP, isotopeIdentification="CAMERA", 
                                    subdir=NULL, nSlaves=2, ppm=15, maxcharge=2)
  checkTrue(all(resultPPCamera$best_settings$result[1:4]== c(0, 221, 65, 101)))
  checkEqualsNumeric(resultPPCamera$best_settings$result[5], 156.9385, tolerance=1e-3)
  
  
  #checking single parameter peak picking optimization
  paramsPP$max_peakwidth <- 15
  paramsPP$ppm <- 50
  resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir=NULL, nSlaves=2, 
                              checkPeakShape="borderIntensity") 
  
  checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 216, 93, 61)))
  checkEqualsNumeric(resultPP$best_settings$result[5], 40.0175, tolerance=1e-3)  
  
  #checking retention time correction and grouping optimization
  mtbls2files <- list.files(paste(find.package("mtbls2"), "/mzData", sep=""), 
                            full.names=TRUE)
  xset <- xcmsSet(mtbls2files[1:2], method="centWave", 
                  peakwidth=c(12, 30), ppm=30, noise=10000, nSlaves=2)
  
  
  #checking obiwarp
  paramsRG <- getDefaultRetGroupStartingParams()
  paramsRG$gapInit <- 0.34
  paramsRG$minfrac <- 1
  paramsRG$profStep <- 1
  paramsRG$mzwid <- 0.026
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir=NULL, nSlaves=2)

  TV <- unlist(resultRG[[2]]$target_value)
  checkTrue(all(TV[-c(5)]== c(1, 185, 0, 34596, 1)))
  checkEqualsNumeric(TV[5], 1723.209, tolerance=1e-2)

  
  #checking single parameter retention time correction and grouping optimization
  paramsRG$gapExtend <- 2.4
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir=NULL, nSlaves=2)
  
  TV <- unlist(resultRG[[2]]$target_value)
  checkTrue(all(TV[-c(5)]== c(1, 185, 0, 34596, 1)))
  checkEqualsNumeric(TV[5], 1723.209, tolerance=1e-2)
  
  #checking loess
  paramsRG <- getDefaultRetGroupStartingParams("loess")
  paramsRG$extra <- 0
  paramsRG$missing <- 0
  resultRG <- optimizeRetGroup(xset, paramsRG, subdir=NULL, nSlaves=2)
  
  TV <- unlist(resultRG[[2]]$target_value)
  checkTrue(all(TV[-c(5)]== c(1, 185, 0, 34596, 1)))
  checkEqualsNumeric(TV[5], 2380, tolerance=1.5e-2)
  
}
