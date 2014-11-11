resultIncreased <-
function(history) {

  index = length(history)
  if(index < 2)
    return(TRUE)
   
  if(history[[index-1]]$max_settings[1] >= history[[index]]$max_settings[1])
    return(FALSE)
    
  return(TRUE)

}
