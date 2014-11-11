createModel <-
function(design, params, resp) {
  design$resp <- resp
  formula <- as.formula(paste("resp ~ SO(", paste("x", 1:length(params), sep="", collapse=","), ")", sep="")) 
  return(rsm(formula, data=design)) 
}
