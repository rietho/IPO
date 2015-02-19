\name{findIsotopes.IPO}
\alias{findIsotopes.IPO}
\title{
Calculation of a peak picking score (PPS) by using natural, stable 13C isotopes
}
\description{
	This function identifies natural, stable 13C isotopes within an xcmsSet object of
	LC-HRMS data. Isotopes have to be within a mass-, retentiontime- and intensity-
	window to be recognized as isotopes. 	
}
\usage{
findIsotopes.IPO(xset)
}

\arguments{
  \item{xset}{\code{xcmsSet} object}
}
\details{
  Identification of 13C isotopes
}
\value{
	An matrix with 2 columns. Column one shows the peak id of the 12C, peak column two 
	shows the id of the respective 13C isotope peak.
}
\author{
	Gunnar Libiseller
}
\seealso{
  \code{\link{findIsotopes.CAMERA}}
}
\examples{

\dontrun{
mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
xset <- xcmsSet(mzmlfile, peakwidth=c(5,12), method="centWave")
isotopes <- findIsotopes.IPO(xset)
}
}

\keyword{method}
