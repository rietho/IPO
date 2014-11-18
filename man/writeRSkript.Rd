\name{writeRSkript}
\alias{writeRSkript}
\title{
	Prints an R-script to the screen which can be used for xcms processing
}
\description{
	This function that prints a script out of the optimized findPeaks, 
	retcor and grouping parameters to the screen.
}
\usage{
writeRSkript(peakPickingSettings, retCorGroupSettings, nSlaves)
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

}
\value{
	none
}
\references{
}
\author{
Gunnar Libiseller
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

	setwd("/...") #set working directory to where files are located\cr
	peakPickingOptimization <- optimizeXcmsSet()
	retcorGroupOptimization <- optimizeRetGroup(peakPickingOptimization$best_settings$xset)
    writeRSkript(peakPickingOptimization$best_settings$parameters, retcorGroupOptimization$best_settings, 4)

## The function is currently defined as
function (peakPickingSettings, retCorGroupSettings, nSlaves) 
{
    cat("library(xcms)\n")
    cat("library(Rmpi)\n")
    if (is.null(peakPickingSettings$step)) {
        cat(paste("xset <- xcmsSet(method=\"centWave\", peakwidth=c(", 
            peakPickingSettings$min_peakwidth, ", ", peakPickingSettings$max_peakwidth, 
            "), ppm=", peakPickingSettings$ppm, ", noise=", peakPickingSettings$noise, 
            ", snthresh=", peakPickingSettings$snthresh, ", mzdiff=", 
            peakPickingSettings$mzdiff, ", prefilter=c(", peakPickingSettings$prefilter, 
            ", ", peakPickingSettings$value_of_prefilter, "), mzCenterFun=\"", 
            peakPickingSettings$mzCenterFun, "\", integrate=", 
            peakPickingSettings$integrate, ", fitgauss=", peakPickingSettings$fitgauss, 
            ", verbose.columns=", peakPickingSettings$verbose.columns, 
            ", nSlaves=", nSlaves, ")\n", sep = ""))
    }
    else {
        cat(paste("xset <- xcmsSet(method=\"matchedFilter\", fwhm=", 
            peakPickingSettings$fwhm, ", snthresh=", peakPickingSettings$snthresh, 
            ", step=", peakPickingSettings$step, ", steps=", 
            peakPickingSettings$steps, ", sigma=", peakPickingSettings$sigma, 
            ", max=", peakPickingSettings$max, ", mzdiff=", peakPickingSettings$mzdiff, 
            ", index=", peakPickingSettings$index, ")\n", sep = ""))
    }
    cat(paste("xset <- retcor(xset, method=\"obiwarp\", plottype=\"", 
        retCorGroupSettings$plottype, "\", distFunc=\"", retCorGroupSettings$distFunc, 
        "\", profStep=", retCorGroupSettings$profStep, ", center=", 
        retCorGroupSettings$center, ", response=", retCorGroupSettings$response, 
        ", gapInit=", retCorGroupSettings$gapInit, ", gapExtend=", 
        retCorGroupSettings$gapExtend, ", factorDiag=", retCorGroupSettings$factorDiag, 
        ", factorGap=", retCorGroupSettings$factorGap, ", localAlignment=", 
        retCorGroupSettings$localAlignment, ")\n", sep = ""))
    cat(paste("xset <- group(xset, method=\"density\", bw=", 
        retCorGroupSettings$bw, ", mzwid=", retCorGroupSettings$mzwid, 
        ", minfrac=", retCorGroupSettings$minfrac, ", minsamp=", 
        retCorGroupSettings$minsamp, ", max=", retCorGroupSettings$max, 
        ")\n", sep = ""))
    cat(paste("xset <- fillPeaks(xset, nSlaves=", nSlaves, ")\n", 
        sep = ""))
  }
}

\keyword{method}
