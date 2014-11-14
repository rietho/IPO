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
}
\value{
	A Central-Composite Design of Experiments
}
\references{
}
\author{
Gunnar Libiseller
}
\note{
}
\seealso{
 \code{\link{getBbdParameter}}
}
\examples{

params <- getDefaultXcmsSetStartingParams()
typ_params <- typeCastParams(params)
design <- getCcdParameter(typ_params$to_optimize) 

## The function is currently defined as
function (params) 
{
    require(rsm)
    lower_bounds <- unlist(lapply(X = params, FUN = function(x) x[1]))
    higher_bounds <- unlist(lapply(X = params, FUN = function(x) x[2]))
    steps <- (higher_bounds - lower_bounds)/2
    x <- paste("x", 1:length(params), " ~ (", c(names(params)), 
        " - ", (lower_bounds + steps), ")/", steps, sep = "")
    formulae <- list()
    for (i in 1:length(x)) formulae[[i]] <- as.formula(x[i])
    design <- ccd(length(params), n0 = 1, alpha = "face", randomize = FALSE, 
        inscribed = TRUE, coding = formulae)
    return(design)
  }
}
\keyword{method}
