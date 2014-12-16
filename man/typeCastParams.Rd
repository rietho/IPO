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
}
\value{
	A list of:
	\item{to_optimize}{A LIST containing all parameters which should be optimized.} 
    \item{no_optimization}{A LIST containing all parameters which should not be optimized.} 
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

## The function is currently defined as
function (params) 
{
    ret_1 <- list()
    ret_2 <- list()
    ret <- list()
    for (i in 1:length(params)) {
        factor <- params[[i]]
        if (length(factor) == 2) {
            ret_1[[(length(ret_1) + 1)]] <- factor
            names(ret_1)[length(ret_1)] <- names(params)[i]
        }
        else {
            ret_2[[(length(ret_2) + 1)]] <- factor
            names(ret_2)[length(ret_2)] <- names(params)[i]
        }
    }
    ret$to_optimize <- ret_1
    ret$no_optimization <- ret_2
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{method}
