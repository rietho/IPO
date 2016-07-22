\name{typeCastParams}
\alias{typeCastParams}
\title{
	Splits parameters into those which should be optimized and those which are fixed.
}
\description{
	This method takes a list of parameters and returns a list consisting of another two lists;
	one holding parameters ment for optimization and one holding fixed parameters.
}
\usage{
	typeCastParams(params)
}

\arguments{
  \item{params}{
	A list of parameters for an xcms-method
}
}
\details{
	This method takes a list of parameters and returns a list consisting of another two lists;
	one holding parameters ment for optimization and one holding fixed parameters.
}
\value{
	A list of:
	\item{to_optimize}{A LIST containing all parameters which should be optimized.} 
    \item{no_optimization}{A LIST containing all parameters which should not be optimized.} 
}
\author{
Gunnar Libiseller
}
\seealso{
\code{\link{optimizeXcmsSet}}, \code{\link{optimizeRetGroup}}
}
\examples{

  params <- getDefaultXcmsSetStartingParams()
  typ_params <- typeCastParams(params)

}
\keyword{method}
