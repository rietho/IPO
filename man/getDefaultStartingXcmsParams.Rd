\name{getDefaultStartingXcmsParams}
\alias{getDefaultStartingXcmsParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getDefaultStartingXcmsParams(method = "centWave")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (method = "centWave") 
{
    if (method == "centWave") 
        return(list(min_peakwidth = c(10, 30), max_peakwidth = c(35, 
            65), ppm = c(15, 35), mzdiff = c(-0.001, 0.01), snthresh = 10, 
            noise = 0, prefilter = 3, value_of_prefilter = 100, 
            mzCenterFun = "wMean", integrate = 1, fitgauss = FALSE, 
            verbose.columns = FALSE))
    if (method == "matchedFilter") 
        return(list(fwhm = c(25, 35), snthresh = c(3, 17), step = c(0.05, 
            0.15), steps = c(1, 3), sigma = 0, max = 5, mzdiff = 0, 
            index = FALSE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
