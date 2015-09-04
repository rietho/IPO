\name{findIsotopes.IPO}
\alias{findIsotopes.IPO}
\title{
Identification of 13C isotopes
}
\description{
	This function identifies natural, stable 13C isotopes within an xcmsSet 
	object of LC-HRMS data. Isotopes have to be within a mass-, retentiontime- 
	and intensitywindow to be recognized as isotopes. If checkBorderIntensity 
	is TRUE the maximum intensity of each peaks has to be at least three times 
	the intensity at rtmin and rtmax. 
}
\usage{
findIsotopes.IPO(xset, checkPeakShape=c("none", "borderIntensity", 
                                        "sinusCurve", "normalDistr"))
}

\arguments{
  \item{xset}{\code{xcmsSet} object}
  \item{checkPeakShape}{
    character to choose if the peakshape should be checked and if so how
  }
}
\details{
  Identification of 13C isotopes
}
\value{
	An matrix with 2 columns. Column one shows the peak id of the 12C, peak 
	column two shows the id of the respective 13C isotope peak.
}
\author{
	Gunnar Libiseller
}
\seealso{
  \code{\link{findIsotopes.CAMERA}}
}
\examples{

mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
xset <- xcmsSet(mzmlfile, peakwidth=c(5,12), method="centWave")
isotopes <- findIsotopes.IPO(xset, "borderIntensity")

}

\keyword{method}
