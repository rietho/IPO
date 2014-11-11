RCSandGSIncreased <-
function(history) {

  index = length(history)
  if(index < 2)
    return(TRUE)
    
  prev_tv <- history[[index-1]]$target_value
  cur_tv <- history[[index]]$target_value
  
  if(cur_tv$bad_groups == 0) {
    cur_tv$bad_groups = 1
    cur_tv$good_groups = cur_tv$good_groups + 1
  }
  
  if(prev_tv$bad_groups == 0) {
    prev_tv$bad_groups = 1
    prev_tv$good_groups = prev_tv$good_groups + 1
  }
  
  if((cur_tv$good_groups^2/cur_tv$bad_groups <= prev_tv$good_groups^2/prev_tv$bad_groups) | (cur_tv$mean_rel_rt_diff >= prev_tv$mean_rel_rt_diff))
    return(FALSE)
    
  return(TRUE)

}
