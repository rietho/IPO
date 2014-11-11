checkRetGroupSetParams <-
function(params) {

  quantitative_parameters <- c("profStep", "gapInit", "gapExtend", "response", "factorDiag", "factorGap", "minfrac", "minsamp", "bw", "mzwid", "max")
  qualitative_parameters <- c("distFunc", "plottype", "localAlignment", "center")
  unsupported_parameters <- c("col", "ty", "initPenalty", "sleep")

  checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
}
