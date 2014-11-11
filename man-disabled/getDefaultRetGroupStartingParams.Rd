\name{getDefaultRetGroupStartingParams}
\alias{getDefaultRetGroupStartingParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getDefaultRetGroupStartingParams(distfunc = "cor_opt", high_resolution = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distfunc}{
%%     ~~Describe \code{distfunc} here~~
}
  \item{high_resolution}{
%%     ~~Describe \code{high_resolution} here~~
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
function (distfunc = "cor_opt", high_resolution = TRUE) 
{
    ret <- NULL
    if (!is.null(distfunc)) {
        if (distfunc == "cor") 
            ret <- (list(distFunc = "cor", gapInit = c(0, 0.4), 
                gapExtend = c(2.1, 2.7)))
        if (distfunc == "cor_opt") 
            ret <- (list(distFunc = "cor_opt", gapInit = c(0, 
                0.4), gapExtend = c(2.1, 2.7)))
        if (distfunc == "cov") 
            ret <- (list(distFunc = "cov", gapInit = c(0, 0.4), 
                gapExtend = c(11.4, 12)))
        if (distfunc == "prd") 
            ret <- (list(distFunc = "prd", gapInit = c(0, 0.4), 
                gapExtend = c(7.5, 8.1)))
        if (distfunc == "euc") 
            ret <- (list(distFunc = "euc", gapInit = c(0.7, 1.1), 
                gapExtend = c(1.5, 2.1)))
        ret$profStep <- c(0.7, 1)
        ret$plottype <- "none"
        ret$response <- 1
        ret$factorDiag <- 2
        ret$factorGap <- 1
        ret$localAlignment <- 0
    }
    else {
        ret <- list()
    }
    ret$bw <- c(22, 38)
    ret$minfrac <- c(0.3, 0.7)
    ret$mzwid <- c(0.15, 0.35)
    if (high_resolution) 
        ret$mzwid <- c(0.015, 0.035)
    ret$minsamp <- 1
    ret$max <- 50
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
