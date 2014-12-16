\name{combineParams}
\alias{combineParams}
\title{
	Combining two lists of parameters together.
}
\description{
	This function combines two lists of parameters. The first is a list of parameters
	which should be optimized. This parameters have different values set by Design of 
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

## The function is currently defined as
function (params_1, params_2) 
{
    len <- max(unlist(sapply(params_1, length)))
    p_names <- c(names(params_1), names(params_2))
    matchedFilter <- "fwhm" \%in\% p_names
    for (i in 1:length(params_2)) {
        new_index <- length(params_1) + 1
        fact <- params_2[[i]]
        params_1[[new_index]] <- fact
        if (matchedFilter) {
            if (p_names[new_index] == "sigma" && fact == 0) {
                if ("fwhm" \%in\% names(params_1)) {
                  params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
                }
                else {
                  params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
                }
            }
            else {
                if (p_names[new_index] == "mzdiff" && fact == 
                  0) {
                  if ("step" \%in\% names(params_1)) {
                    if ("steps" \%in\% names(params_1)) {
                      params_1[[new_index]][1:len] <- 0.8 - params_1$step * 
                        params_1$steps
                    }
                    else {
                      params_1[[new_index]][1:len] <- 0.8 - params_1$step * 
                        params_2$steps
                    }
                  }
                  else {
                    if ("steps" \%in\% names(params_1)) {
                      params_1[[new_index]][1:len] <- 0.8 - params_2$step * 
                        params_1$steps
                    }
                    else {
                      params_1[[new_index]][1:len] <- 0.8 - params_2$step * 
                        params_2$steps
                    }
                  }
                }
                else {
                  params_1[[new_index]][1:len] <- fact
                }
            }
        }
        else {
            params_1[[new_index]][1:len] <- fact
        }
    }
    names(params_1) <- p_names
    return(params_1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{method}
