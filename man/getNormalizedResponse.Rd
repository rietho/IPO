\name{getNormalizedResponse}
\alias{getNormalizedResponse}
\title{
	It combines Retention time Correction Scores (RCS) and Grouping Scores (GS)
}
\description{
	This function does unity based normalization on Retention time Correction 
	Scores (RCS) as well as Grouping Scores (GS). 
}
\usage{
getNormalizedResponse(response)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
	A List of all responses calculated by getRGTVValues for all experiments of an
	Design of Experiment
}
}
\details{
	Grouping Score (GS)	is calculated by:\cr
	'good groups'^2/'bad groups\cr\cr
	
	For all RCS and GS values unitiy based normalization is done. For every 
	experiment within the DoE these two values are added together and returned.
}
\value{
	A vector with RTGV values 
}
\author{
Gunnar Libiseller
}
\note{
	Since RCS and GS can be within completely different ranges, normalization has
	to be done to prevent an excessive influence of either RCS or GS.
}

\seealso{
  \code{\link{getRGTVValues}}
}
\examples{
\donttest{

mtbls2files <- list.files(file.path(find.package("mtbls2"), "mzML"), 
                          full.names=TRUE)

params <- list(min_peakwidth=12, max_peakwidth=30, ppm=30,
                  mzdiff=-0.001, snthresh=10, noise=10000, prefilter=3, 
                  value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
                  fitgauss=FALSE, verbose.columns=FALSE, nSlaves=2)

xset <- calculateXcmsSet(mtbls2files[1:2], params)
xset <- retcor(xset, method="obiwarp")
xset <- group(xset)

result <- getRGTVValues(xset)
result
}
}

\keyword{method}
