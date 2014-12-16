\name{createModel}
\alias{createModel}
\title{
	Creating a response surface model.
}
\description{
	This function uses a design of experiments, a response for the experiments within the
	design and the used parameters to create a response surface model
}
\usage{
createModel(design, params, resp)
}

\arguments{
  \item{design}{
	A design of experiments (Box-Behnken-Design or Central-Composite-Design)
}
  \item{params}{
	The parameters which were used.
}
  \item{resp}{
	The responses achivied for the various experiments.
}
}
\details{
	This function uses a design of experiments, a response for the experiments within the
	design and the used parameters to create a response surface model
}
\value{
	A response surface model.
}
\references{
Lenth, R. V. (2009). Response-Surface Methods in R , Using rsm. Journal
  of Statistical Software, 32(7), 1-17. Retrieved from \url{http://www.jstatsoft.org/v32/i07}
}
\author{
	Gunnar Libiseller
}
\note{
  \code{\link{getBbdParameter}}
  \code{\link{getCcdParameter}}
  \code{\link{typeCastParams}}
}

\examples{

  params <- getDefaultXcmsSetStartingParams()
  type_params <- typeCastParams(params)
  design <- getBbdParameter(type_params$to_optimize)
  resp <- runif(nrow(design),1,3)
  
  model <- createModel(design, type_params$to_optimize, resp)
  dev.new()
  par(mfrow=c(3,2))
  contour(model, ~ x1*x2*x3*x4, image=TRUE)

  

}

\keyword{method}