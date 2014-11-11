\name{typeCastFactor}
\alias{typeCastFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
typeCastFactor(params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
