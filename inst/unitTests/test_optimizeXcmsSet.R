test_ipo <- function() {
	mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
 
	paramsPP <- getDefaultXcmsSetStartingParams()
	paramsPP$mzdiff <- -0.001
	paramsPP$min_peakwidth <- c(7,14)
	paramsPP$max_peakwidth <- c(20,30)
	resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir="microtofq")

    	checkTrue(all(resultPP$best_settings$result[1:4]== c(0, 246, 83, 96)))
        checkEqualsNumeric(resultPP$best_settings$result[5], 111.0361, tolerance=0.000001)
}
