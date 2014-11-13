\name{getNormalizedResponse}
\alias{getNormalizedResponse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNormalizedResponse(response)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
%%     ~~Describe \code{response} here~~
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
function (response) 
{
    good_groups <- sapply(response, "[[", "good_groups")
    bad_groups <- sapply(response, "[[", "bad_groups")
    bad_groups_bool <- bad_groups == 0
    bad_groups[bad_groups_bool] <- 1
    good_groups[bad_groups_bool] <- good_groups[bad_groups_bool] + 
        1
    group_ratio <- good_groups^2/bad_groups
    ARTS <- 1/sapply(response, "[[", "mean_rel_rt_diff")
    ARTS_penalty <- 1/sapply(response, "[[", "retcor_done")
    ARTS <- ARTS/ARTS_penalty
    norm_group_ratio <- (group_ratio - min(group_ratio))/(max(group_ratio) - 
        min(group_ratio))
    norm_ARTS <- (ARTS - min(ARTS))/(max(ARTS) - min(ARTS))
    return(norm_group_ratio + norm_ARTS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
