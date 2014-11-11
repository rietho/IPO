getCcdParameter <-
function(params) {
 
  require(rsm)
  lower_bounds <- unlist(lapply(X=params, FUN=function(x) x[1]))
  higher_bounds <- unlist(lapply(X=params, FUN=function(x) x[2]))
  
  steps <- (higher_bounds - lower_bounds)/2
  
  x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", (lower_bounds + steps), ")/", steps, sep="")
  formulae <- list()
  for(i in 1:length(x))
    formulae[[i]] <- as.formula(x[i])  
  
  design <- ccd(length(params), n0 = 1, alpha = "face", randomize = FALSE, inscribed = TRUE, coding = formulae)

  return(design)
  
}
