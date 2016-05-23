\name{writeParamsTable}
\alias{writeParamsTable}
\title{
	Writes XCMS settings to a file.
}
\description{
	This function writes findPeaks, retcor and grouping parameters to a
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
	This function writes findPeaks, retcor and grouping parameters to a
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

\examples{

#creating list of peak picking parameters
paramsPP <- list(min_peakwidth=5, max_peakwidth=12, ppm=58,
                  mzdiff=-0.001, snthresh=10, noise=0, prefilter=3, 
                  value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
                  fitgauss=FALSE, verbose.columns=FALSE, nSlaves=1)
                  

#creating list of retention time correction and grouping parameters
paramsRTCGroup <- list(retcorMethod="obiwarp", distFunc="cor", gapInit=0.2, 
                       gapExtend=2.4, profStep=1, plottype="none", response=1,  
                       factorDiag=2, factorGap=1, localAlignment=0, initPenalty=0,
                       bw=30, minfrac=0.5, minsamp=1, mzwid=0.25, max=50)

#writing parameters to the file "params.tsv"
writeParamsTable(paramsPP, paramsRTCGroup, "params.tsv")


}

\keyword{method}
