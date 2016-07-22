\name{IPO-package}
\alias{IPO-package}
\alias{IPO}
\docType{package}
\title{Automated Optimization of Untargeted Metabolomics LC-MS Data Processing}
\description{
IPO provides a framework for parameter optimization for the software package XCMS.
It provides optimisation of peak picking parameters by using natural, stable 13C 
isotopes. Retention time correction is optimized by minimizing the relative retention
time differences within features and grouping parameters are optimized by maximizing
the number of features showing exactly one peak from each injection of a pooled sample.
}
\details{
An overview of how to use the package, including the most important functions
}
\author{
Gunnar Libiseller

Maintainer: gunnar.libiseller@joanneum.at
}
\references{

  Lenth, R. V. (2009). Response-Surface Methods in R , Using rsm. Journal of Statistical Software, 32(7), 1-17. Retrieved from \url{http://www.jstatsoft.org/v32/i07}

  Smith, C.A. and Want, E.J. and O'Maille, G. and Abagyan,R. and
  Siuzdak, G.: XCMS: Processing mass spectrometry data for metabolite
  profiling using nonlinear peak alignment, matching and
  identification, Analytical Chemistry, 78:779-787 (2006)

  Ralf Tautenhahn, Christoph Boettcher, Steffen Neumann: Highly
  sensitive feature detection for high resolution LC/MS BMC
  Bioinformatics, 9:504 (2008)

  H. Paul Benton, Elizabeth J. Want and Timothy M. D. Ebbels Correction
  of mass calibration gaps in liquid chromatography-mass spectrometry
  metabolomics data Bioinformatics, 26:2488 (2010)
  
  Yu, H. (2002). Rmpi: Parallel Statistical Computing in R. R News,
  2(2), 10-14. 
  Retrieved from \url{http://cran.r-project.org/doc/Rnews/Rnews_2002-2.pdf}

}
\seealso{
\code{\link[xcms:xcmsSet]{xcms}}
}
\examples{
\dontrun{
mtbls2files <- list.files(file.path(find.package("mtbls2"), "mzData"), 
                          full.names=TRUE)

 
paramsPP <- getDefaultXcmsSetStartingParams()
paramsPP$mzdiff <- -0.001
#paramsPP$ppm <- 25
paramsPP$min_peakwidth <- c(7,14)
paramsPP$max_peakwidth <- c(20,30)
paramsPP$noise <- 10000
resultPP <- optimizeXcmsSet(mtbls2files[1:2], paramsPP, subdir="mtbls2")


paramsRG <- getDefaultRetGroupStartingParams()
paramsRG$gapInit <- 0.2
paramsRG$profStep <- 1
paramsRG$minfrac <- 0.75
resultRG <- optimizeRetGroup(resultPP$best_settings$xset, paramsRG, nSlaves=2)

writeRScript(resultPP$best_settings$parameters, resultRG$best_settings, 
             subdir="mtbls2", 4)
}
}
