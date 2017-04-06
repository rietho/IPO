\name{calculateXcmsSet}
\alias{calculateXcmsSet}
\title{
Calculation of an xcmsSet-Object
}
\description{
	This function encapsulates xcms::findPeaks-methods for IPO
}
\usage{
calculateXcmsSet(files, xcmsSetParameters, scanrange=NULL, task=1, 
  BPPARAM = bpparam(), nSlaves=0) 
}

\arguments{
  \item{files}{a vector containing the files for peak picking}
  \item{xcmsSetParameters}{
    a list with all parameters for \code{\link[xcms]{xcmsSet}}-methods as list-items}
  \item{scanrange}{
     scan range to process. See \code{\link[xcms]{findPeaks.centWave}}.}
  \item{task}{
    The task-id when using this method in parallel calculations.}
  \item{BPPARAM}{a \code{BiocParallel} parameter object to control how
    and if parallel processing of \code{\link[xcms]{xcmsSet}} should be performed. 
    Such objects can be created by the \code{\link[BiocParallel]{SerialParam}},
    \code{\link[BiocParallel]{MulticoreParam}} or
    \code{\link[BiocParallel]{SnowParam}} functions.
  }
  \item{nSlaves}{
    \code{\link[xcms]{xcmsSet}}'s \code{nSlaves}-argument is deprecated., use \code{BPPARAM} argument instead.}
}
\details{
  Encapsulation of xcms::findPeaks-methods used in IPO.
}
\value{
	An xcmsSet-object
}
\references{

  Smith, C.A. and Want, E.J. and O'Maille, G. and Abagyan,R. and
  Siuzdak, G.: XCMS: Processing mass spectrometry data for metabolite
  profiling using nonlinear peak alignment, matching and
  identification, Analytical Chemistry, 78:779-787 (2006)

  Ralf Tautenhahn, Christoph Boettcher, Steffen Neumann: Highly
  sensitive feature detection for high resolution LC/MS BMC
  Bioinformatics, 9:504 (2008)
}
\author{
	Gunnar Libiseller,
	Thomas Riebenbauer (\email{thomas.riebenbauer@joanneum.at})
}
\seealso{
\code{\link[xcms]{findPeaks}}
}
\examples{

mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")

params <- list(min_peakwidth=5, max_peakwidth=12, ppm=58,
                  mzdiff=-0.001, snthresh=10, noise=0, prefilter=3, 
                  value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
                  fitgauss=FALSE, verbose.columns=FALSE, nSlaves=1)

xset <- calculateXcmsSet(mzmlfile, params)

}

\keyword{method}
