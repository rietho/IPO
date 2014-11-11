\name{optimizeRetGroupSlave}
\alias{optimizeRetGroupSlave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
optimizeRetGroupSlave()
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
    while (done != 1) {
        mpi.send.Robj(junk, 0, 1)
        task <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
        task_info <- mpi.get.sourcetag()
        tag <- task_info[2]
        print(paste("tag", tag, "task", task))
        if (tag == 1) {
            require(xcms)
            exp_index <- task
            do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && 
                is.null(parameters$gapInit) && is.null(parameters$gapExtend) && 
                is.null(parameters$plottype) && is.null(parameters$response) && 
                is.null(parameters$factorDiag) && is.null(parameters$factorGap) && 
                is.null(parameters$localAlignment) && is.null(parameters$initPenalty))
            retcor_failed = ifelse(do_retcor, 1.1, 1)
            if (do_retcor) {
                try(retcor_failed <- retcor(xset, method = "obiwarp", 
                  plottype = parameters$plottype[exp_index], 
                  distFunc = parameters$distFunc[exp_index], 
                  profStep = parameters$profStep[exp_index], 
                  center = parameters$center[exp_index], response = parameters$response[exp_index], 
                  gapInit = parameters$gapInit[exp_index], gapExtend = parameters$gapExtend[exp_index], 
                  factorDiag = parameters$factorDiag[exp_index], 
                  factorGap = parameters$factorGap[exp_index], 
                  localAlignment = parameters$localAlignment[exp_index]))
                if (!is.numeric(retcor_failed)) {
                  xset <- retcor_failed
                  retcor_failed = 1
                }
            }
            minfrac <- ifelse(is.null(parameters$minfrac), 1, 
                parameters$minfrac[exp_index])
            try(xset <- group(xset, method = "density", bw = parameters$bw[exp_index], 
                mzwid = parameters$mzwid[exp_index], minfrac = minfrac, 
                minsamp = parameters$minsamp[exp_index], max = parameters$max[exp_index]))
            result <- calculateRGTV(xset, exp_index, retcor_failed)
            mpi.send.Robj(result, 0, 2)
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
