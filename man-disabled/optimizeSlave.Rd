\name{optimizeSlave}
\alias{optimizeSlave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optimizeSlave()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    junk <- 0
    done <- 0
    library(Rmpi)
    while (done != 1) {
        mpi.send.Robj(junk, 0, 1)
        task <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
        task_info <- mpi.get.sourcetag()
        tag <- task_info[2]
        if (tag == 1) {
            require(xcms)
            exp_index <- task
            if (is.null(xcmsSet_parameters$step)) {
                print(sapply(xcmsSet_parameters, "[[", exp_index))
                xset <- xcmsSet(files = example_sample, method = "centWave", 
                  peakwidth = c(xcmsSet_parameters$min_peakwidth[exp_index], 
                    xcmsSet_parameters$max_peakwidth[exp_index]), 
                  ppm = xcmsSet_parameters$ppm[exp_index], noise = xcmsSet_parameters$noise[exp_index], 
                  snthresh = xcmsSet_parameters$snthresh[exp_index], 
                  mzdiff = xcmsSet_parameters$mzdiff[exp_index], 
                  prefilter = c(xcmsSet_parameters$prefilter[exp_index], 
                    xcmsSet_parameters$value_of_prefilter[exp_index]), 
                  mzCenterFun = xcmsSet_parameters$mzCenterFun[exp_index], 
                  integrate = xcmsSet_parameters$integrate[exp_index], 
                  fitgauss = xcmsSet_parameters$fitgauss[exp_index], 
                  verbose.columns = xcmsSet_parameters$verbose.columns[exp_index])
            }
            else {
                xset <- xcmsSet(files = example_sample, method = "matchedFilter", 
                  fwhm = xcmsSet_parameters$fwhm[exp_index], 
                  snthresh = xcmsSet_parameters$snthresh[exp_index], 
                  step = xcmsSet_parameters$step[exp_index], 
                  steps = xcmsSet_parameters$steps[exp_index], 
                  sigma = xcmsSet_parameters$sigma[exp_index], 
                  max = xcmsSet_parameters$max[exp_index], mzdiff = xcmsSet_parameters$mzdiff[exp_index], 
                  index = xcmsSet_parameters$index[exp_index])
            }
            result <- calcPPS(xset)
            result[1] <- exp_index
            rm(xset)
            mpi.send.Robj(result, 0, 2)
            print("result sent")
        }
        else if (tag == 2) {
            done <- 1
        }
    }
    mpi.send.Robj(junk, 0, 3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
