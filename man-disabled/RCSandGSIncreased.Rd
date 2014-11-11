\name{RCSandGSIncreased}
\alias{RCSandGSIncreased}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RCSandGSIncreased(history)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{history}{
%%     ~~Describe \code{history} here~~
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
function (history) 
{
    index = length(history)
    if (index < 2) 
        return(TRUE)
    prev_tv <- history[[index - 1]]$target_value
    cur_tv <- history[[index]]$target_value
    if (cur_tv$bad_groups == 0) {
        cur_tv$bad_groups = 1
        cur_tv$good_groups = cur_tv$good_groups + 1
    }
    if (prev_tv$bad_groups == 0) {
        prev_tv$bad_groups = 1
        prev_tv$good_groups = prev_tv$good_groups + 1
    }
    if ((cur_tv$good_groups^2/cur_tv$bad_groups <= prev_tv$good_groups^2/prev_tv$bad_groups) | 
        (cur_tv$mean_rel_rt_diff >= prev_tv$mean_rel_rt_diff)) 
        return(FALSE)
    return(TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
