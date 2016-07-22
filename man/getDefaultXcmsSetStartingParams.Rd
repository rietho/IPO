\name{getDefaultXcmsSetStartingParams}
\alias{getDefaultXcmsSetStartingParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Creates a List of parameters for xcms-methods xcmsSet.findPeak which are optimized by default
}
\description{
  This function creates a list of parameters used in the xcmsSet.findPeak-methods 
  'centWave' and 'matchedFilter'. Per default the following parameters 
  have a defined range where optimization should start:\cr
  'centWave' parameters: 'peakwidth' (split into 'min_peakwidth' and 'max_peakwidth'), 
  'ppm', 'mzdiff'\cr
  'matchedFilter' parameters: 'fwhm', 'snthresh', 'step', 'steps'\cr  
}
\usage{
getDefaultXcmsSetStartingParams(method = c("centWave", "matchedFilter"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{Either parameters for 'centWave' or 'matchedFilter' should be created
}
}
\details{
  * Do not delete a parameter from the list returned.\cr
  * Optimization of qualitative parameters is not supported yet.\cr
  * If you want to optimize additional parameter just set an lower and
  an upper bound (e.g. params$snthresh <- c(5,20))\cr
  * If you dont want to optimize a parameter set a default value 
  (e.g. params$snthresh <- 10)
}
\value{
  A List of parameters for the xcmsSet.findPeak-methods 
  'centWave' or 'matchedFilter'
}
\author{
Gunnar Libiseller
}
\examples{
params <- getDefaultXcmsSetStartingParams()
params$ppm <- 10
params$snthresh <- c(5,15)
params

params <- getDefaultXcmsSetStartingParams("matchedFilter")
params

}
\keyword{method}

