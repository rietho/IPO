\name{decode}
\alias{decode}
\title{
	Decodes a values from ranges of -1 to 1 to ranges specified.
}
\description{
	A function used to decode values that are in a range of -1 to 1 into a  
	specified range.
}
\usage{
	decode(value, bounds)
}

\arguments{
  \item{value}{
	A value 
}
  \item{bounds}{
	A vector of two values defining the lower and upper bound of a range.
}
}
\details{
}
\value{
	The encoded value.
}
\references{
}
\author{
	Gunnar Libiseller
}
\note{
}

\seealso{
  \code{\link{decodeAll}}
  \code{\link{encode}}
}
\examples{

	decode(0, c(10, 20))
	decode(-0.5, c(10, 20))
	decode(1, c(10, 20))
	
	bounds <- c(10, 20)
	encode(decode(1, bounds), bounds)

}
\keyword{method}
