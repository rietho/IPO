\name{writeRSkript}
\alias{writeRSkript}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeRSkript(peakPickingSettings, retCorGroupSettings, nSlaves)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peakPickingSettings}{
%%     ~~Describe \code{peakPickingSettings} here~~
}
  \item{retCorGroupSettings}{
%%     ~~Describe \code{retCorGroupSettings} here~~
}
  \item{nSlaves}{
%%     ~~Describe \code{nSlaves} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
