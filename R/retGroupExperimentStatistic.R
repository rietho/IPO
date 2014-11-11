retGroupExperimentStatistic <-
function(retcor_result, subdir, iterator, xset) {

  params <- retcor_result$params
  resp <- getNormalizedResponse(retcor_result$response)
  model <- createModel(retcor_result$design, params$to_optimize, resp)
  
  retcor_result$model <- model                  
  max_settings <- getMaximumExperiment(retcor_result$model)
  tmp <- max_settings[-1]
  tmp[is.na(tmp)] <- 1

  plotContours(retcor_result$model, tmp, paste(subdir, "/retgroup_rsm_", iterator, sep=""))  
    
  parameters <- as.list(decodeAll(max_settings[-1], params$to_optimize)) 
  parameters <- combineParams(parameters, params$no_optimization)
  xset_tmp <- xset
  
  exp_index <- 1				   
  do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && is.null(parameters$gapInit) && is.null(parameters$gapExtend)
	             && is.null(parameters$plottype) && is.null(parameters$col) && is.null(parameters$ty) && is.null(parameters$response) 
				 && is.null(parameters$factorDiag) && is.null(parameters$factorGap) && is.null(parameters$localAlignment) 
				 && is.null(parameters$initPenalty)) #&& is.null(parameters$center)) 
					
      retcor_failed = ifelse(do_retcor, 1.1, 1)  
  
  if(do_retcor) {
    try(retcor_failed <- retcor(xset_tmp, method="obiwarp", plottype=parameters$plottype[exp_index], distFunc=parameters$distFunc[exp_index],
                             profStep=parameters$profStep[exp_index], center=parameters$center[exp_index],
							 response=parameters$response[exp_index], gapInit=parameters$gapInit[exp_index], 
				 	         gapExtend=parameters$gapExtend[exp_index], factorDiag=parameters$factorDiag[exp_index], 
					         factorGap=parameters$factorGap[exp_index], localAlignment=parameters$localAlignment[exp_index]))
  	
	
    if(!is.numeric(retcor_failed)) {
      xset_tmp <- retcor_failed
      retcor_failed=1
    } 
  }
	  
  minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])      
  try(xset_tmp <- group(xset_tmp, method="density", bw=parameters$bw[exp_index], mzwid=parameters$mzwid[exp_index], minfrac=minfrac, 
                  minsamp=parameters$minsamp[exp_index], max=parameters$max[exp_index]))	 
 
  tv <- calculateRGTV(xset_tmp, exp_index, retcor_failed)

  retcor_result$max_settings <- max_settings
  retcor_result$target_value <- tv   
  retcor_result$xset <- xset_tmp
  retcor_result$best_settings <- parameters

  return(retcor_result)

}
