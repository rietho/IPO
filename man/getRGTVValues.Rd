\name{calculateRGTV}
\alias{calculateRGTV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateRGTV(xset, exp_index = 1, retcor_penalty = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
  \item{exp_index}{
%%     ~~Describe \code{exp_index} here~~
}
  \item{retcor_penalty}{
%%     ~~Describe \code{retcor_penalty} here~~
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
function (xset, exp_index = 1, retcor_penalty = 1) 
{
    features <- nrow(xset@groups)
    good_groups <- sum(unlist(lapply(X = xset@groupidx, FUN = function(x, 
        xset) {
        ifelse(length(unique(xset@peaks[x, "sample"])) == length(xset@filepaths) & 
            length(xset@peaks[x, "sample"]) == length(xset@filepaths), 
            1, 0)
    }, xset)))
    bad_groups <- nrow(xset@groups) - good_groups
    relative_rt_diff <- c()
    if (nrow(xset@groups) > 0) {
        for (i in 1:nrow(xset@groups)) {
            feature_rtmed <- xset@groups[i, "rtmed"]
            relative_rt_diff <- c(relative_rt_diff, mean(abs(feature_rtmed - 
                xset@peaks[xset@groupidx[[i]], "rt"])/feature_rtmed))
        }
    }
    else {
        relative_rt_diff <- 1
    }
    ARTS <- (mean(relative_rt_diff)) * retcor_penalty
    ret <- list(exp_index = exp_index, good_groups = good_groups, 
        bad_groups = bad_groups, mean_rel_rt_diff = ARTS)
    ret$retcor_done = retcor_penalty
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
