\name{checkParams}
\alias{checkParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{quantitative_parameters}{
%%     ~~Describe \code{quantitative_parameters} here~~
}
  \item{qualitative_parameters}{
%%     ~~Describe \code{qualitative_parameters} here~~
}
  \item{unsupported_parameters}{
%%     ~~Describe \code{unsupported_parameters} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (params, quantitative_parameters, qualitative_parameters, 
    unsupported_parameters) 
{
    for (i in 1:length(params)) {
        param <- params[[i]]
        name <- names(params)[i]
        if (name \%in\% unsupported_parameters) {
            stop(paste("The parameter", name, "is not supported! Please remove from parameters; stopping!"))
        }
        if (name \%in\% qualitative_parameters) {
            if (length(param) == 0) {
                stop(paste("The parameter", name, "has no value set! Please specify; stopping!"))
            }
            if (length(param) > 1) {
                stop(paste("Optimization of parameter", name, 
                  "no supported! Please specify only one value; stopping!"))
            }
        }
        if (name \%in\% quantitative_parameters) {
            if (length(param) == 0) {
                stop(paste("The parameter", name, "has no value set! Please specify between one and two; stopping!"))
            }
            if (length(param) > 2) {
                stop(paste("Too many values for parameter", name, 
                  "! Please specify only one or two; stopping!"))
            }
        }
    }
    missing_params <- which(!(c(quantitative_parameters, qualitative_parameters) \%in\% 
        names(params)))
    if (length(missing_params > 0)) {
        stop(paste("The parameter(s)", paste(c(quantitative_parameters, 
            qualitative_parameters)[missing_params], collapse = ", "), 
            "are missing! Please specify; stopping!"))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
