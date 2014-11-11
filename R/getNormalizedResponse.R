getNormalizedResponse <-
function(response) {

  good_groups <- sapply(response, "[[", "good_groups")
  bad_groups <- sapply(response, "[[", "bad_groups")
  bad_groups_bool <- bad_groups == 0
  bad_groups[bad_groups_bool] <- 1
  good_groups[bad_groups_bool] <- good_groups[bad_groups_bool] + 1
  group_ratio <- good_groups ^ 2 / bad_groups 
  ARTS <- 1/sapply(response, "[[", "mean_rel_rt_diff")
  
  #give penalty when retcor failed
  ARTS_penalty <- 1/sapply(response, "[[", "retcor_done")
  ARTS <- ARTS/ARTS_penalty
  
  #normalize
  norm_group_ratio <- (group_ratio - min(group_ratio)) / (max(group_ratio) - min(group_ratio))
  norm_ARTS <- (ARTS - min(ARTS)) / (max(ARTS) - min(ARTS))
  
  return(norm_group_ratio+norm_ARTS)

}
