\name{getCcdParameter}
\alias{getCcdParameter}
\title{
	Creates a Central-Composite Design of experiment
}
\description{
	Creates a Central-Composite Design of Experiment out of a list of parameters.
	Each of the list items has to be a pair defining the lower und upper limits
	of the value-range to test. The method then returns a Center faced Central-Composite
	Design of Experiments. 
	
}
\usage{
getCcdParameter(params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
	A list of value pairs defining lower und upper limits of an optimization
	range.
}
}
\details{
	Creates a Central-Composite Design of Experiment out of a list of parameters.
	Each of the list items has to be a pair defining the lower und upper limits
	of the value-range to test. The method then returns a Center faced Central-Composite
	Design of Experiments. 
}
\value{
	A Central-Composite Design of Experiments
}
\references{
  Lenth, R. V. (2009). Response-Surface Methods in R , Using rsm. Journal of Statistical Software, 32(7), 1-17. Retrieved from 
  \url{http://www.jstatsoft.org/v32/i07}
}
\author{
Gunnar Libiseller
}
\seealso{
 \code{\link{getBbdParameter}}
}
\examples{
params <- getDefaultXcmsSetStartingParams()
typ_params <- typeCastParams(params)
design <- getCcdParameter(typ_params$to_optimize) 
}
\keyword{method}
