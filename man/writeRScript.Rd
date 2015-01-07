\name{writeRScript}
\alias{writeRScript}
\title{
	Prints an R-script to the screen which can be used for xcms processing
}
\description{
	This function prints a script of the optimized findPeaks, 
	retcor and grouping parameters to the screen.
}
\usage{
writeRScript(peakPickingSettings, retCorGroupSettings, nSlaves)
}

\arguments{
  \item{peakPickingSettings}{
	The optimized settings for xcms-methods findPeaks.centWave or findPeaks.matchedFilter
}
  \item{retCorGroupSettings}{
	The optimized settings for xcms-methods for retcor.obiwarp and group.density
}
  \item{nSlaves}{
	The number of slaves which should be used in the script
}
}
\details{
	This function that prints a script out of the optimized findPeaks, 
	retcor and grouping parameters to the screen.
  
}
\value{
	none
}
\author{
Gunnar Libiseller
}

\keyword{method}
