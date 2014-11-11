\name{plotContours}
\alias{plotContours}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotContours(model, maximum_slice, plot_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{maximum_slice}{
%%     ~~Describe \code{maximum_slice} here~~
}
  \item{plot_name}{
%%     ~~Describe \code{plot_name} here~~
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
function (model, maximum_slice, plot_name) 
{
    plots <- c()
    for (i in 1:(length(maximum_slice) - 1)) {
        for (j in (i + 1):length(maximum_slice)) {
            plots <- c(plots, as.formula(paste("~ x", i, "* x", 
                j, sep = "")))
        }
    }
    plot_rows <- round(sqrt(length(plots)))
    plot_cols <- if (plot_rows == 1) {
        length(plots)
    }
    else {
        ceiling(sqrt(length(plots)))
    }
    if (!is.null(plot_name)) {
        plot_name = paste(plot_name, ".jpg", sep = "")
        jpeg(plot_name, width = 4 * plot_cols, height = 2 * plot_rows + 
            2, unit = "in", res = c(200, 200))
    }
    else dev.new(width = 4 * plot_cols, height = 2 * plot_rows + 
        2)
    par(mfrow = c(plot_rows, plot_cols), oma = c(3, 0, 2, 0))
    contours <- contour(model, plots, image = TRUE, at = maximum_slice)
    if (!is.null(plot_name)) {
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
