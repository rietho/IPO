\name{combineParams}
\alias{combineParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combineParams(params_1, params_2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params_1}{
%%     ~~Describe \code{params_1} here~~
}
  \item{params_2}{
%%     ~~Describe \code{params_2} here~~
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
function (params_1, params_2) 
{
    len <- max(unlist(sapply(params_1, length)))
    p_names <- c(names(params_1), names(params_2))
    matchedFilter <- "fwhm" \%in\% p_names
    for (i in 1:length(params_2)) {
        new_index <- length(params_1) + 1
        fact <- params_2[[i]]
        params_1[[new_index]] <- fact
        if (matchedFilter) {
            if (p_names[new_index] == "sigma" && fact == 0) {
                if ("fwhm" \%in\% names(params_1)) {
                  params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
                }
                else {
                  params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
                }
            }
            else {
                if (p_names[new_index] == "mzdiff" && fact == 
                  0) {
                  if ("step" \%in\% names(params_1)) {
                    if ("steps" \%in\% names(params_1)) {
                      params_1[[new_index]][1:len] <- 0.8 - params_1$step * 
                        params_1$steps
                    }
                    else {
                      params_1[[new_index]][1:len] <- 0.8 - params_1$step * 
                        params_2$steps
                    }
                  }
                  else {
                    if ("steps" \%in\% names(params_1)) {
                      params_1[[new_index]][1:len] <- 0.8 - params_2$step * 
                        params_1$steps
                    }
                    else {
                      params_1[[new_index]][1:len] <- 0.8 - params_2$step * 
                        params_2$steps
                    }
                  }
                }
                else {
                  params_1[[new_index]][1:len] <- fact
                }
            }
        }
        else {
            params_1[[new_index]][1:len] <- fact
        }
    }
    names(params_1) <- p_names
    return(params_1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
