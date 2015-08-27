\name{calcPPS}
\alias{calcPPS}
\title{
Calculation of a peak picking score (PPS) by using natural, stable 13C isotopes
}
\description{
	This function calculates PPS by identifying natural, stable 13C isotopes within
	an xcmsSet object. Peaks beeing part of an isotopologue are defined
	as reliable peaks (RP). Peaks which are not part of an isotopologue
	and where the intensity of possible isotopes is below a cutoff are defined as 
	'low intensity peaks' (LIP). PPS is then calculated by:
	\cr
	PPS = RP^2/(#all_peaks - LIP)
}
\usage{
calcPPS(xset, isotopeIdentification, ...)
}

\arguments{
  \item{xset}{\code{xcmsSet} object}
  \item{isotopeIdentification}{This parameter defines the method for isotope 
        identification. The method IPO was especially implemented for high resolution data. 
        CAMERA is an established isotope and adduct annotation package.}
  \item{...}{
     Additional parameters to CAMERA's findIsotopes function.}
}
\details{
  Calculation of a peak picking score (PPS) by using natural, stable 13C isotopes
}
\value{
	An array with 5 items:\cr
	1. Space for experimentid of the Design of Experiments (0 since not known in calcPPS)\cr
	2. Number of peaks\cr
	3. Number of peaks without LIP and RP\cr
	4. Reliable peaks (RP)\cr
	5. Peak picking score (PPS)\cr
}
\author{
	Gunnar Libiseller
}
\seealso{
  \code{\link{findIsotopes.IPO}}
  \code{\link{findIsotopes.CAMERA}}
}
\examples{
mzmlfile <- file.path(find.package("msdata"), "microtofq/MM14.mzML")
xset <- xcmsSet(mzmlfile, peakwidth=c(5,12), method="centWave")
calcPPS(xset)
}

\keyword{method}
