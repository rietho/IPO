retGroupCalcExperiments <-
function(params, xset, n_slaves=4) {
							   
  library(Rmpi)  
  library(rsm)
    
  junk <- 0
  closed_slaves <- 0
  #n_slaves <- min(mpi.comm.size()-1, n_slaves)  
  
  typ_params <- typeCastFactor(params)
  
  if(length(typ_params$to_optimize) > 2) {
    design <- getBbdParameter(typ_params$to_optimize) 
  } else {
    design <- getCcdParameter(typ_params$to_optimize) 
  }	
 
  parameters <- decode.data(design)	
  tasks <- as.list(1:nrow(design))    
  startSlaves(n_slaves)
  
  parameters <- combineParams(parameters, typ_params$no_optimization)
  
  sendRetGroupSlaveFunctions(parameters, xset) 
  
  response <- list()  
  finished <- 0
  while(closed_slaves < n_slaves) {
    # Receive a message from a slave
    message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    message_info <- mpi.get.sourcetag()
    slave_id <- message_info[1]
    tag <- message_info[2]
      
    if(tag == 1) {
      if(length(tasks) > 0) {
        mpi.send.Robj(tasks[[1]], slave_id, 1);
        tasks[[1]] <- NULL
      } else {
        mpi.send.Robj(junk, slave_id, 2)
      }
    } else if (tag == 2) {
	  #print(message)
      response[[message$exp_index]] <- message
      finished <- finished + 1     
    } else if (tag == 3) {
      # A slave has closed down. 
      closed_slaves <- closed_slaves + 1
    }
    cat(paste("finished ", finished, " of ", nrow(design), " tasks\r", sep="")) 
    flush.console()
  }
  cat("\n\r")  
  print("done")

  mpi.close.Rslaves()
  
  ret <- list()
  ret$params <- typ_params
  ret$design <- design
  #ret$model <- model
  ret$response <- response
  
  return(ret)

}
