\name{checkXcmsSetParams}
\alias{checkXcmsSetParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkXcmsSetParams(params)
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
    if (is.null(params$step)) {
        quantitative_parameters <- c("ppm", "min_peakwidth", 
            "max_peakwidth", "snthresh", "mzdiff", "noise", "prefilter", 
            "value_of_prefilter")
        qualitative_parameters <- c("integrate", "fitgauss", 
            "verbose.columns", "mzCenterFun")
        unsupported_parameters <- c("scanrange", "sleep", "ROI.list")
    }
    else {
        quantitative_parameters <- c("fwhm", "sigma", "max", 
            "snthresh", "step", "steps", "mzdiff")
        qualitative_parameters <- c("index")
        unsupported_parameters <- c("sleep")
    }
    checkParams(params, quantitative_parameters, qualitative_parameters, 
        unsupported_parameters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
