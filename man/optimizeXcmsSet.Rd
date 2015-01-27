\name{optimizeXcmsSet}
\alias{optimizeXcmsSet}
\title{
Optimisation of peak picking parameters by using natural, stable 13C isotopes
}
\description{
    This function provides optimisation of peak picking parameters by using natural, stable 13C isotopes.
}
\usage{
    optimizeXcmsSet(files, params = getDefaultXcmsSetStartingParams(),
nSlaves = 4, subdir = "IPO")
}
\arguments{
  \item{files}{A directory or list of files, passed to xcmsSet(). If no
    files are given, xcmsSet() will check recursively all MS files in
    the current working directory. }
  \item{params}{
     A list of parameters which are needed by XCMS::findPeaks-Methods. List-items with
     two values will be optimized. The first value defines the lower test value, the 
     second one the higher test value.
  }
  \item{nSlaves}{Number of slaves the optimization process should spawn.}
  \item{subdir}{
     The name of the subdirectory which is created and where the figures of the response
     surface models will be saved to.
  }
}
\details{
    This function provides optimisation of peak picking parameters by using natural, stable 13C isotopes.

}
\value{
    A LIST of length n+1 with n beeing the optimization runs needed

  \item{comp1-comp(n) }{A LIST containing:\cr 
    1. Parameters used for the nth optimization run\cr 
	2. Box-Behnken or Central Composite Design used for optimization run\cr 
	3. Responses from calcPPS for every experiment in the design\cr 
	4. Response surface model for the design\cr 
	5. The normalized parameter settings giving the best PPS (values between -1 and 1)}
  \item{comp(n+1) }{A LIST containing:\cr   
    1. Parameters giving the best PPS\cr 
	2. An xcmsSet-object created with the optimized parameters\cr 
	3. The result of calcPPS() of the xcmsSet created with the best parameters}
}
\references{

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

}
\author{
Gunnar Libiseller
}
\seealso{
  \code{\link{getDefaultXcmsSetStartingParams}}
  \code{\link{calcPPS}}
}
\examples{
\dontrun{
library(IPO)
 
mzdatafiles <- list.files(file.path(find.package("mtbls2"), "mzdata"), full.names=TRUE)[1:4]
 
paramsPP <- getDefaultXcmsSetStartingParams()
paramsPP$mzdiff <- -0.001
paramsPP$min_peakwidth <- c(7,14)
paramsPP$max_peakwidth <- c(20,30)
resultPP <- optimizeXcmsSet(mzdatafiles, paramsPP, subdir="mtbls2")
}
}
\keyword{method}

