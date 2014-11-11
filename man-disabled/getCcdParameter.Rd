\name{getCcdParameter}
\alias{getCcdParameter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCcdParameter(params)
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
    require(rsm)
    lower_bounds <- unlist(lapply(X = params, FUN = function(x) x[1]))
    higher_bounds <- unlist(lapply(X = params, FUN = function(x) x[2]))
    steps <- (higher_bounds - lower_bounds)/2
    x <- paste("x", 1:length(params), " ~ (", c(names(params)), 
        " - ", (lower_bounds + steps), ")/", steps, sep = "")
    formulae <- list()
    for (i in 1:length(x)) formulae[[i]] <- as.formula(x[i])
    design <- ccd(length(params), n0 = 1, alpha = "face", randomize = FALSE, 
        inscribed = TRUE, coding = formulae)
    return(design)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
