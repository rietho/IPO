\name{getDefaultRetCorCenterSample}
\alias{getDefaultRetCorCenterSample}
	
\title{
	Gets the index of the sample with most peaks in it.
}
\description{
	Gets the index of the sample with most peaks in it. This is used if no center sample
	for retention time correction has been defined by the user.
}
\usage{
getDefaultRetCorCenterSample(xset)
}
\arguments{
  \item{xset}{\code{xcmsSet} object}
}
\details{
	Gets the index of the sample with most peaks in it. This is used if no center sample
	for retention time correction has been defined by the user.
}
\value{
	The file index of the sample with most peaks in it.
}
\author{
Gunnar Libiseller
}
\examples{


## The function is currently defined as
function (xset) 
{
    ret <- NULL
    for (i in 1:length(filepaths(xset))) {
        ret <- c(ret, sum(peaks(xset)[, "sample"] == i))
    }
    return(which.max(ret))
  }
}

\keyword{method}
