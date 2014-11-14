\name{attachList}
\alias{attachList}
\title{
	Attaching one list at the end of another
}
\description{
	This function attaches one list at the end of another list.
}
\usage{
params1 <- attachList(params_1, params_2)
}
\arguments{
  \item{params_1}{
	A List
}
  \item{params_2}{
	A second list which will be attached at the end of the first list.
}
}
\details{

}
\value{
	A List composed of the two input lists.
}
\references{

}
\author{
Gunnar Libiseller
}
\note{

}

\seealso{
}
\examples{
	a <- list("a"=1, "b"=2)
	b <- list("c"=4, "d"=4)
	attachList(a, b)
}
\keyword{method}
