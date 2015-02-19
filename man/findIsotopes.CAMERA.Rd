\name{findIsotopes.CAMERA}
\alias{findIsotopes.CAMERA}
\title{
Calculation of a peak picking score (PPS) by using natural, stable 13C isotopes
}
\description{
	This function finds isotopes using CAMERA's find peak function. Isotopes are 
	separately found within each sample. 	
}
\usage{
findIsotopes.CAMERA(xset, ...)
}

\arguments{
  \item{xset}{\code{xcmsSet} object}
  \item{...}{Additional parameters to the findIsotopes function of CAMERA}
}
\details{
  Identification of 13C isotopes
}
\value{
	An matrix with 2 columns. Column one shows the peak id of the 12C, peak column two 
	shows the id of the respective 13C isotope peak. 
}
\references{
  C. Kuhl and R. Tautenhahn and C. Boettcher and T. R. Larson and S. Neumann: CAMERA: 
  an integrated strategy for compound spectra extraction
  and annotation of liquid chromatography/mass spectrometry data sets 
  Analytical Chemistry 84:283 (2012)
}
\author{
	Gunnar Libiseller
}
\seealso{
  \code{\link{findIsotopes.IPO}}
}
\examples{

\dontrun{
mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
xset <- xcmsSet(mzmlfile, peakwidth=c(5,12), method="centWave")
isotopes <- findIsotopes.CAMERA(xset, ppm=15, maxcharge=2)
}
}

\keyword{method}
