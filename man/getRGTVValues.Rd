\name{getRGTVValues}
\alias{getRGTVValues}
\title{
	Calculation of values used for a Retention time correction and Grouping Target Value (RGTV)
}
\description{
	This function calculates the Average Retention Time Shifts (ARTS) of all features within
	an xcmsSet-object. Also features having exactly one peak from each sample are defined as
	'good groups', all others a 'bad groups'.
}
\usage{
getRGTVValues(xset, exp_index = 1, retcor_penalty = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset}{
{\code{xcmsSet} object}
}
  \item{exp_index}{
	Experiment-id of the experiment within a Design of Experiments
}
  \item{retcor_penalty}{
	Penalty if an error occured with the used retention time correction parameters
}
}
\details{
	This function calculates the Average Retention Time Shifts (ARTS) of all features within
	an xcmsSet-object. Also features having exactly one peak from each sample are defined as
	'good groups', all others a 'bad groups'.
}
\value{
  something.
}
\author{
Gunnar Libiseller
}
\examples{
\dontrun{
getRGTVValues(xset)
}
}
\keyword{method}
