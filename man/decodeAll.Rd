\name{decodeAll}
\alias{decodeAll}
\title{
	Decodes a vector of values from ranges of -1 to 1 to ranges specified.
}
\description{
	A function used to decode values that are in a range of -1 to 1 into a  
	specified range. For every value a list item with lower and upper bound
	has to be supplied.
}
\usage{
    decodeAll(values, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{
	A vector with values in the range [-1,1]
}
  \item{params}{
	A list where evere list-item consist of two values defining a lower
	and an upper bound.
}
}
\details{
}
\value{
	A vector of decoded values.
}
\references{
}
\author{
  Gunnar Libiseller
}
\note{
}

\seealso{
  \code{\link{decode}}
  \code{\link{encode}}
}
\examples{
	values <- c(-1, -0.25, 0, 0.75)
	params <- getDefaultXcmsSetStartingParams()
    type_params <- typeCastParams(params)
	
	decodeAll(values, type_params$to_optimize)
}

\keyword{methdo}
