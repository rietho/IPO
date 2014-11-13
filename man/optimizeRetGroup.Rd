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
optimizeRetGroup(xset, params = getDefaultRetGroupStartingParams(), n_slaves = 4, subdir = "IPO")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset}{
	\code{xcmsSet} object
}
  \item{params}{
     A list of parameters which are needed by xcms-methods retcor.obiwarp and group.density. 
	 List-items with two values will be optimized. The first value defines the lower test value, the 
     second one the higher test value.
  }
  \item{n_slaves}{
Number of slaves the optimization process should spawn.
  }
  \item{subdir}{
The name of the subdirectory which is created and where the figures of the response
     surface models will be saved to.
  }
}
\details{
}
\value{
    A LIST of length n+1 with n beeing the optimization runs needed

  \item{comp1-comp(n) }{A LIST containing:\cr 
    1. Parameters used for the nth optimization run\cr 
	2. Box-Behnken or Central Composite Design used for optimization run\cr 
	3. Responses from calculateRGTV for every experiment in the design\cr 
	4. Response surface model for the design\cr 
	5. The normalized parameter settings giving the best Retention time and Grouping Target Value (RGTV) (values between -1 and 1)}
	6. Response from calculateRGTV for xcmsSet-object created with best parameters in this run
	7. xcmsSet-object created with best parameters in this run
  \item{comp(n+1) }{A LIST containing:\cr   
    1. Parameters giving the best RGTV\cr 

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gunnar Libiseller
}
\note{
}
\seealso{
  \code{\link{getDefaultRetGroupStartingParams}
  \code{\link{calculateRGTV)}
}
\examples{
	setwd("/...") #set working directory to where files are located\cr
	peakPickingOptimization <- optimizeXcmsSet()
	retcorGroupOptimization <- optimizeRetGroup(peakPickingOptimization$best_settings$xset)
}
\keyword{method}
