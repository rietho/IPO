resultIncreased <-
function(history) {

  index = length(history)
  if(history[[index]]$max_settings[1] == 0 & index == 1)
    stop("No isotopes have been detected, peak picking not optimizable by IPO!")
  
  if(index < 2)
    return(TRUE)
   
  if(history[[index-1]]$max_settings[1] >= history[[index]]$max_settings[1])
    return(FALSE)
    
  return(TRUE)

}
