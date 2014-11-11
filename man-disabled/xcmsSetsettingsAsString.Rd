\name{xcmsSetsettingsAsString}
\alias{xcmsSetsettingsAsString}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xcmsSetsettingsAsString(parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
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
function (parameters) 
{
    return(sprintf("xcmsSet(method='centWave', ppm=\%s, peakwidth=c(\%s,\%s), snthresh=\%s, prefilter=c(\%s,\%s), mzCenterFun='\%s', integrate=\%s, mzdiff=\%s, fitgauss=\%s, noise=\%s, verbose.columns=\%s)", 
        prettyNum(parameters$ppm), prettyNum(parameters$min_peakwidth), 
        prettyNum(parameters$max_peakwidth), prettyNum(parameters$snthresh), 
        prettyNum(parameters$prefilter), prettyNum(parameters$value_of_prefilter), 
        parameters$mzCenterFun, prettyNum(parameters$integrate), 
        prettyNum(parameters$mzdiff), parameters$fitgauss, prettyNum(parameters$noise), 
        parameters$verbose.columns))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
