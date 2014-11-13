\name{optimizeXcmsSet}
\alias{optimizeXcmsSet}
\title{
Optimisation of peak picking parameters by using natural, stable 13C isotopes
}
\description{
    This function provides optimisation of peak picking parameters by using natural, stable 13C isotopes.
}
\usage{
    optimizeXcmsSet(params = getDefaultStartingXcmsParams(), n_slaves = 4, subdir = "IPO")
}
\arguments{
  \item{params}{
     A list of parameters which are needed by XCMS::findPeaks-Methods. List-items with
     two values will be optimized. The first value defines the lower test value, the 
     second one the higher test value.
  }
  \item{n_slaves}{Number of slaves the optimization process should spawn.}
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
	3. Responses from calcPPS for every experiment in the design\cr 
	4. Response surface model for the design\cr 
	5. The normalized parameter settings giving the best PPS (values between -1 and 1)}
  \item{comp(n+1) }{A LIST containing:\cr   
    1. Parameters giving the best PPS\cr 
	2. An xcmsSet-object created with the optimized parameters\cr 
	3. The result of calcPPS() of the xcmsSet created with the best parameters}
}
\references{
% ~put references to the literature/web site here ~
}
\author{
Gunnar Libiseller
}
\note{

}
\seealso{
  \code{\link{getDefaultStartingXcmsParams}
  \code{\link{calcPPS)}
}
}
\examples{
	setwd("/...") #set working directory to where files are located\cr
	peakPickingOptimization <- optimizeXcmsSet()
}
\keyword{method}

