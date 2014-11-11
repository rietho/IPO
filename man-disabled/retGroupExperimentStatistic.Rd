\name{retGroupExperimentStatistic}
\alias{retGroupExperimentStatistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
retGroupExperimentStatistic(retcor_result, subdir, iterator, xset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{retcor_result}{
%%     ~~Describe \code{retcor_result} here~~
}
  \item{subdir}{
%%     ~~Describe \code{subdir} here~~
}
  \item{iterator}{
%%     ~~Describe \code{iterator} here~~
}
  \item{xset}{
%%     ~~Describe \code{xset} here~~
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
function (retcor_result, subdir, iterator, xset) 
{
    params <- retcor_result$params
    resp <- getNormalizedResponse(retcor_result$response)
    model <- createModel(retcor_result$design, params$to_optimize, 
        resp)
    retcor_result$model <- model
    max_settings <- getMaximumExperiment(retcor_result$model)
    tmp <- max_settings[-1]
    tmp[is.na(tmp)] <- 1
    plotContours(retcor_result$model, tmp, paste(subdir, "/retgroup_rsm_", 
        iterator, sep = ""))
    parameters <- as.list(decodeAll(max_settings[-1], params$to_optimize))
    parameters <- combineParams(parameters, params$no_optimization)
    xset_tmp <- xset
    exp_index <- 1
    do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && 
        is.null(parameters$gapInit) && is.null(parameters$gapExtend) && 
        is.null(parameters$plottype) && is.null(parameters$col) && 
        is.null(parameters$ty) && is.null(parameters$response) && 
        is.null(parameters$factorDiag) && is.null(parameters$factorGap) && 
        is.null(parameters$localAlignment) && is.null(parameters$initPenalty))
    retcor_failed = ifelse(do_retcor, 1.1, 1)
    if (do_retcor) {
        try(retcor_failed <- retcor(xset_tmp, method = "obiwarp", 
            plottype = parameters$plottype[exp_index], distFunc = parameters$distFunc[exp_index], 
            profStep = parameters$profStep[exp_index], center = parameters$center[exp_index], 
            response = parameters$response[exp_index], gapInit = parameters$gapInit[exp_index], 
            gapExtend = parameters$gapExtend[exp_index], factorDiag = parameters$factorDiag[exp_index], 
            factorGap = parameters$factorGap[exp_index], localAlignment = parameters$localAlignment[exp_index]))
        if (!is.numeric(retcor_failed)) {
            xset_tmp <- retcor_failed
            retcor_failed = 1
        }
    }
    minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])
    try(xset_tmp <- group(xset_tmp, method = "density", bw = parameters$bw[exp_index], 
        mzwid = parameters$mzwid[exp_index], minfrac = minfrac, 
        minsamp = parameters$minsamp[exp_index], max = parameters$max[exp_index]))
    tv <- calculateRGTV(xset_tmp, exp_index, retcor_failed)
    retcor_result$max_settings <- max_settings
    retcor_result$target_value <- tv
    retcor_result$xset <- xset_tmp
    retcor_result$best_settings <- parameters
    return(retcor_result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
