\name{retGroupCalcExperiments}
\alias{retGroupCalcExperiments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
retGroupCalcExperiments(params, xset, n_slaves = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
  \item{n_slaves}{
%%     ~~Describe \code{n_slaves} here~~
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
function (params, xset, n_slaves = 4) 
{
    library(Rmpi)
    library(rsm)
    junk <- 0
    closed_slaves <- 0
    typ_params <- typeCastFactor(params)
    if (length(typ_params$to_optimize) > 2) {
        design <- getBbdParameter(typ_params$to_optimize)
    }
    else {
        design <- getCcdParameter(typ_params$to_optimize)
    }
    parameters <- decode.data(design)
    tasks <- as.list(1:nrow(design))
    startSlaves(n_slaves)
    parameters <- combineParams(parameters, typ_params$no_optimization)
    sendRetGroupSlaveFunctions(parameters, xset)
    response <- list()
    finished <- 0
    while (closed_slaves < n_slaves) {
        message <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
        message_info <- mpi.get.sourcetag()
        slave_id <- message_info[1]
        tag <- message_info[2]
        if (tag == 1) {
            if (length(tasks) > 0) {
                mpi.send.Robj(tasks[[1]], slave_id, 1)
                tasks[[1]] <- NULL
            }
            else {
                mpi.send.Robj(junk, slave_id, 2)
            }
        }
        else if (tag == 2) {
            response[[message$exp_index]] <- message
            finished <- finished + 1
        }
        else if (tag == 3) {
            closed_slaves <- closed_slaves + 1
        }
        cat(paste("finished ", finished, " of ", nrow(design), 
            " tasks\r", sep = ""))
        flush.console()
    }
    cat("\n\r")
    print("done")
    mpi.close.Rslaves()
    ret <- list()
    ret$params <- typ_params
    ret$design <- design
    ret$response <- response
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
