\name{optimizeRetGroup}
\alias{optimizeRetGroup}
\title{
	Optimization for parameters for retention time correction and grouping
}
\description{
	This function provides optimisation for parameters of the xcms-method retcor.obiwarp and
	group.density. The retention time correction is optimised by minimizing intra-feature 
	retention time shifts; grouping is optimized by increasing the number of features which
	have exactly one peak per sample.
}
\usage{
  optimizeRetGroup(xset, params = getDefaultRetGroupStartingParams(),
                   nSlaves = 4, subdir = "IPO")
}

\arguments{
  \item{xset}{
	\code{xcmsSet} object
}
  \item{params}{
     A list of parameters which are needed by xcms-methods retcor.obiwarp and 
     group.density. List-items with two values will be optimized. The first value 
     defines the lower test value, the second one the higher test value.
  }
  \item{nSlaves}{
Number of slaves the optimization process should spawn.
  }
  \item{subdir}{
The name of the subdirectory which is created and where the figures of the response
     surface models will be saved to. \code{NULL} plots the figures to the graphic device.
  }
}
\details{
	This function provides optimisation for parameters of the xcms-method retcor.obiwarp and
	group.density. The retention time correction is optimised by minimizing intra-feature 
	retention time shifts; grouping is optimized by increasing the number of features which
	have exactly one peak per sample.
}
\value{
    A LIST of length n+1 with n beeing the optimization runs needed

  \item{comp1-comp(n) }{A LIST containing:\cr 
  1. Parameters used for the nth optimization run\cr 
	2. Box-Behnken or Central Composite Design used for optimization run\cr 
	3. Responses from calculateRGTV for every experiment in the design\cr 
	4. Response surface model for the design\cr 
	5. The normalized parameter settings giving the best Retention time and 
     Grouping Target Value (RGTV) (values between -1 and 1)
	6. Response from calculateRGTV for xcmsSet-object created with best parameters 
     in this run
	7. xcmsSet-object created with best parameters in this run
  }
  \item{comp(n+1) }{A LIST containing:\cr   
    1. Parameters giving the best RGTV\cr 
  }
}
\references{

Obiwarp
  \code{Prince, J. T., & Marcotte, E. M. (2006). Chromatographic alignment of 
  ESI-LC-MS proteomics data sets by ordered bijective interpolated warping. 
  Analytical chemistry, 78(17), 6140-52. doi:10.1021/ac0605344}

}
\author{
Gunnar Libiseller
}
\seealso{
  \code{\link{getDefaultRetGroupStartingParams}}
%  \code{\link{calculateRGTV}}
}
\examples{

\donttest{
 

mtbls2files <- list.files(file.path(find.package("mtbls2"), "mzData"), 
                          full.names=TRUE)

params <- list(min_peakwidth=12, max_peakwidth=30, ppm=30,
                  mzdiff=-0.001, snthresh=10, noise=10000, prefilter=3, 
                  value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
                  fitgauss=FALSE, verbose.columns=FALSE, nSlaves=2)

xset <- calculateXcmsSet(mtbls2files[1:2], params)


#optimize the retention time correction and grouping parameters
paramsRG <- getDefaultRetGroupStartingParams()
paramsRG$profStep <- 1
paramsRG$minfrac <- 0.75
resultRG <- optimizeRetGroup(xset, params=paramsRG, 
                             nSlaves=4,subdir="mtbls2")

writeRScript(params, resultRG$best_settings, 4)

}
}
\keyword{method}
