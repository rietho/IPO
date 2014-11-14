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
	data <- toMatrix(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	An array or a matrix
}
}
\details{

}
\value{
  A matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gunnar Libsieller
}
\note{
}

\seealso{
}
\examples{
	data <- matrix(1:9, nrow=3)\cr
	colnames(data) <- c("a","b","c")\cr
	x <- data[1,]
	is.matrix(x)
	x <- toMatrix(x)
	is.matrix(x)
}
\keyword{method}
