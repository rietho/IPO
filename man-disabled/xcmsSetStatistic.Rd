\name{xcmsSetStatistic}
\alias{xcmsSetStatistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xcmsSetStatistic(xcms_result, subdir, iterator, score_name = "PPS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcms_result}{
%%     ~~Describe \code{xcms_result} here~~
}
  \item{subdir}{
%%     ~~Describe \code{subdir} here~~
}
  \item{iterator}{
%%     ~~Describe \code{iterator} here~~
}
  \item{score_name}{
%%     ~~Describe \code{score_name} here~~
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
function (xcms_result, subdir, iterator, score_name = "PPS") 
{
    params <- xcms_result$params
    resp <- xcms_result$response[, "PPS"]
    model <- createModel(xcms_result$design, params$to_optimize, 
        resp)
    xcms_result$model <- model
    max_settings <- getMaximumExperiment(xcms_result$model)
    tmp <- max_settings[-1]
    tmp[is.na(tmp)] <- 1
    plotContours(xcms_result$model, tmp, paste(subdir, "/rsm_", 
        iterator, sep = ""))
    xcms_result$max_settings <- max_settings
    return(xcms_result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
