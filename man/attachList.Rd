\name{attachList}
\alias{attachList}
\title{
	Attaching one list at the end of another
}
\description{
	This function attaches one list at the end of another list.
}
\usage{
attachList(params_1, params_2)
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
  This is a convenience funktion, but the implementation
  is not optimized for speed. 
}
\value{
	A List composed of the two input lists.
}
\author{
Gunnar Libiseller
}
\examples{
	a <- list("a"=1, "b"=2)
	b <- list("c"=4, "d"=4)
	attachList(a, b)
}
\keyword{method}
