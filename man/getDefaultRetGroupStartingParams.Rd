\name{getDefaultRetGroupStartingParams}
\alias{getDefaultRetGroupStartingParams}

\title{
	Gives a List parameters for xcms-methods retcor.obiwarp and group.density 
	which are optimized by default
}
\description{
  This function creates a list of parameters used in the xcms-methods 
  retcor.obiwarp and group.density. Per default the following parameters 
  have a defined range where optimization should start:\cr
  retcor.obiwarp parameters: 'gapInit'; 'gapExtend', 'profStep'\cr
  group.density parameters: 'bw', 'minfrac', 'mzwid'\cr  
}
\usage{
getDefaultRetGroupStartingParams(distfunc = "cor_opt", high_resolution = TRUE)
}

\arguments{
  \item{distfunc}{
	The name of the distance function used by retcor.obiwarp
}
  \item{high_resolution}{
	If high_resolution = TRUE starting values for mzwid are set to 0.015 and 0.035;
	if high_resolution = FALSE to 0.15, 0.35
}
}
\details{
  * Do not delete a parameter from the list returned.\cr
  * Optimization of qualitative parameters is not supported yet.\cr
  * If you want to optimize additional parameter just set an lower and
  an upper bound (e.g. params$max <- c(4,8))\cr
  * If you dont want to optimize a parameter set a default value 
  (e.g. params$max <- 10)
}
\value{
  A List of parameters used in the xcms-methods retcor.obiwarp and group.density
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
params <- getDefaultRetGroupStartingParams()
params$bw <- 10
params$max <- c(4,8)
params
}

\keyword{method}

