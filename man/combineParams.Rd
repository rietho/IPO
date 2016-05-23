\name{combineParams}
\alias{combineParams}
\title{
	Combining two lists of parameters together.
}
\description{
	This function combines two lists of parameters. The first is a list of parameters
	which should be optimized. These parameters have different values set by Design of 
	Experimnt. The second list consists of parameters which should not be optimized, 
	hence only one value is set for each parameter. The parameters of the second list
	are replicated to have the same length as the number experiments in the DoE. Then
	the two lists are combined.
}
\usage{
combineParams(params_1, params_2)
}
\arguments{
  \item{params_1}{
	A list holding parameters which should be optimized. Each parameter already has
	value set for each experiment of an Design of Experiment.
}
  \item{params_2}{	
	A list holding parameters which should not be optimized, hence only one value is set.
}
}
\details{
	Special treatment is needed for the findPeaks.matchedFilter-parameters 'sigma', 'mzdiff'
	since these two parameters are defined by default relative to the parameters 'fwhm' or 
	'step' and 'steps' respectively.\cr
	\code{sigma=fwhm/2.3548}
	\code{mzdiff=0.8-step*steps}
}
\value{
	A list of consting of all parameters needed for an xcms-method (findPeaks.centWave, 
	findPeaks.matchedFilter, retcor.obiwarp or group.density). Each list item has the 
	same length which is equal to the number of experiments within the DoE.
}
\author{
Gunnar Libiseller
}
\examples{

params <- getDefaultXcmsSetStartingParams()
typ_params <- typeCastParams(params)
design <- getBbdParameter(typ_params$to_optimize) 
xcms_design <- decode.data(design) 
xcms_design <- combineParams(xcms_design, typ_params$no_optimization)  
xcms_design

}

\keyword{method}
