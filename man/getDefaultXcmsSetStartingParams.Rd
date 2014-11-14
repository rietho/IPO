\name{getDefaultXcmsSetStartingParams}
\alias{getDefaultXcmsSetStartingParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Gives a List parameters for xcms-methods xcmsSet.findPeak which are optimized by default
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
getDefaultXcmsSetStartingParams(method = "centWave")
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
\references{

}
\author{
Gunnar Libiseller
}
\note{

}



\seealso{
}
\examples{
params <- getDefaultXcmsSetStartingParams()
params$ppm <- 10
params$snthresh <- c(5,15)
params


## The function is currently defined as
function (method = "centWave") 
{
    if (method == "centWave") 
        return(list(min_peakwidth = c(10, 30), max_peakwidth = c(35, 
            65), ppm = c(15, 35), mzdiff = c(-0.001, 0.01), snthresh = 10, 
            noise = 0, prefilter = 3, value_of_prefilter = 100, 
            mzCenterFun = "wMean", integrate = 1, fitgauss = FALSE, 
            verbose.columns = FALSE))
    if (method == "matchedFilter") 
        return(list(fwhm = c(25, 35), snthresh = c(3, 17), step = c(0.05, 
            0.15), steps = c(1, 3), sigma = 0, max = 5, mzdiff = 0, 
            index = FALSE))
  }
}
\keyword{method}

