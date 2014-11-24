\name{encode}
\alias{encode}
\title{
	Encodes a values from range specified into the range of -1 to 1 to .
}
\description{
	A function used to encode values that are in a specified range into a  
	range between -1 to 1.
}
\usage{
	encode(value, bounds)
}

\arguments{
  \item{value}{
	A value in the range [-1, 1]
}
  \item{bounds}{
	A vector of two values defining the lower and upper bound of the new range.
}
}
\details{
}
\value{
	A decoded value.
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
  \code{\link{decodeAll}}
}
\examples{

	encode(15, c(10, 20))
	encode(10, c(10, 20))
	encode(5, c(1, 5))

	bounds <- c(1,5)
	decode(encode(5, bounds), bounds)
}
\keyword{method}
