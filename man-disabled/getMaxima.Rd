\name{getMaxima}
\alias{getMaxima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMaxima(slice, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slice}{
%%     ~~Describe \code{slice} here~~
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
function (slice, model) 
{
    values <- contour(model, ~x1 * x2, at = slice, plot.it = FALSE, 
        decode = FALSE)
    z_values <- values[[1]]$z
    max_z <- max(z_values)
    maxima <- which(z_values == max_z, arr.ind = T)
    ret <- c(max_z, 0, 0)
    if (length(maxima) == 2) {
        ret[2] <- values[[1]]$x[maxima[1]]
        ret[3] <- values[[1]]$y[maxima[2]]
    }
    else {
        for (i in 1:2) {
            col_max <- maxima[, i]
            if (sum(col_max == 1) > 0 && sum(col_max == length(values[[1]][[i]])) > 
                0) 
                ret[(i + 1)] <- NA
            else ret[(i + 1)] <- values[[1]][[i]][max(col_max)]
        }
    }
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
