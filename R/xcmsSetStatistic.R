xcmsSetStatistic <-
function(xcms_result, subdir, iterator, score_name="PPS") {

  params <- xcms_result$params
  resp <- xcms_result$response[, "PPS"]
  model <- createModel(xcms_result$design, params$to_optimize, resp)
  xcms_result$model <- model                  
     
  max_settings <- getMaximumExperiment(xcms_result$model)
  tmp <- max_settings[-1]
  tmp[is.na(tmp)] <- 1
  if(!is.null(subdir))
    plotContours(xcms_result$model, tmp, paste(subdir,"/rsm_", iterator, sep=""))

  xcms_result$max_settings <- max_settings

  return(xcms_result)
}
