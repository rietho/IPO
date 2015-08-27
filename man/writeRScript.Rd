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

#outputting an xcms-script to the display
writeRScript(paramsPP, paramsRTCGroup, 4)


}

\keyword{method}
