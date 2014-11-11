\name{getResponses}
\alias{getResponses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getResponses(slices, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slices}{
%%     ~~Describe \code{slices} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (slices, model) 
{
    add_col <- array(0, dim = c(nrow(slices), 2))
    slices <- cbind(add_col, slices)
    colnames(slices) <- paste("x", 1:ncol(slices), sep = "")
    values <- apply(X = slices, MARGIN = 1, FUN = getMaxima, 
        model)
    ret <- cbind(t(values), slices[, -c(1, 2)])
    colnames(ret) <- c("response", paste("x", 1:ncol(slices), 
        sep = ""))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
