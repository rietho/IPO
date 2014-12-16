\name{toMatrix}
\alias{toMatrix}
\title{
	Converts an array into a matrix
}
\description{
	This function converts an array into a matrix. This is useful to counter the implicit
	casting of matrices into arrays when only one row is selected. If a matrix is passed
	to the function, the matrix is returned, if an array is passed, a matrix with one row
	is returned.
}
\usage{
	toMatrix(data)
}
\arguments{
  \item{data}{
	An array or a matrix
}
}
\value{
  A matrix
}
\author{
Gunnar Libsieller
}
\examples{
	data <- matrix(1:9, nrow=3)
	colnames(data) <- c("a","b","c")
	x <- data[1,]
	is.matrix(x)
	x <- toMatrix(x)
	is.matrix(x)
}
\keyword{method}
