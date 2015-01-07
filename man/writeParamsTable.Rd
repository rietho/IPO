\name{writeParamsTable}
\alias{writeParamsTable}
\title{
	Writes XCMS settings to a file.
}
\description{
	This function write findPeaks, retcor and grouping parameters to a
	file using write.table.
}
\usage{
writeParamsTable(peakPickingSettings, retCorGroupSettings, file, ...)
}

\arguments{
  \item{peakPickingSettings}{
	A list of optimized settings for xcms-methods findPeaks.centWave or findPeaks.matchedFilter
}
  \item{retCorGroupSettings}{
	A list of optimized settings for xcms-methods for retcor.obiwarp and group.density
}
  \item{file}{
	The name of the outputfile for the parameters.
}
  \item{...}{
	Additional parameters for write.table.
}
}
\details{
	This function write findPeaks, retcor and grouping parameters to a
	file using write.table.
  
}
\value{
	none
}
\author{
Gunnar Libiseller
}
\seealso{
\code{\link[xcms:xcmsSet]{xcms}}
}
\keyword{method}
