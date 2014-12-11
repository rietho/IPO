\name{IPO-package}
\alias{IPO-package}
\alias{IPO}
\docType{package}
\title{Automated Optimization of Untargeted Metabolomics LC-MS Data Processing}
\description{
IPO provides a framework for parameter optimization for the software package XCMS.
It provides optimisation of peak picking parameters by using natural, stable 13C isotopes. Retention time correction is optimized by minimizing the relative retention time differences 
within features and grouping parameters are optimized by maximizing the number of features 
showing exactly one peak from each injection of a pooled sample.
}
\details{
An overview of how to use the package, including the most important functions
}
\author{
Gunnar Libiseller

Maintainer: gunnar.libiseller@joanneum.at

}
\references{

  Lenth, R. V. (2009). Response-Surface Methods in R , Using rsm. Journal of Statistical Software, 32(7), 1–17. Retrieved from \link{http://www.jstatsoft.org/v32/i07}

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
  
  Yu, H. (2002). Rmpi: Parallel Statistical Computing in R. R News, 2(2), 10–14. Retrieved from \link{http://cran.r-project.org/doc/Rnews/Rnews_2002-2.pdf}

}
\keyword{
}
\seealso{
\code{\link[xcms:xcmsSet]{xcms}}
}
\examples{
library(msdata)
library(IPO)
 
mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
 
paramsPP <- getDefaultXcmsSetStartingParams()
paramsPP$mzdiff <- -0.001
paramsPP$min_peakwidth <- c(7,14)
paramsPP$max_peakwidth <- c(20,30)
resultPP <- optimizeXcmsSet(mzmlfile, paramsPP, subdir="microtofq")
#resultRG <- optimizeRetGroup(resultPP$best_settings$xset)

#writeRScript(resultPP$best_settings$parameters, resultRG$best_settings, 4)
}
