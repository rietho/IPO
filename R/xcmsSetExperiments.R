xcmsSetExperiments <-
function(example_sample, params, nSlaves=4) { #ppm=5, rt_diff=0.01, nSlaves=4) {
  library(Rmpi)  
  library(rsm)
    
  junk <- 0
  closed_slaves <- 0
  #nSlaves <- min(mpi.comm.size()-1, nSlaves)  
  
  typ_params <- typeCastParams(params)
  
  if(length(typ_params[[1]])>2)
    design <- getBbdParameter(typ_params$to_optimize) 
  else
    design <- getCcdParameter(typ_params$to_optimize)  	
  xcms_design <- decode.data(design) 

  xcms_design <- combineParams(xcms_design, typ_params$no_optimization)  
  tasks <- as.list(1:nrow(design))  
  
  startSlaves(nSlaves)
  sendXcmsSetSlaveFunctions(example_sample, xcms_design) #,  ppm, rt_diff)

  response <- matrix(0, nrow=length(design[[1]]), ncol=5)
  colnames(response) <- c("exp", "num_peaks", "notLLOQP", "num_C13", "PPS")
  finished <- 0
  while(closed_slaves < nSlaves) {
    # Receive a message from a slave
    message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    message_info <- mpi.get.sourcetag()
    slave_id <- message_info[1]
    tag <- message_info[2]
      
    if(tag == 1) {
      # slave is ready for a task.  Give it the next task, or tell it tasks
      # are done if there are none.
      if(length(tasks) > 0) {
        # Send a task, and then remove it from the task list
        mpi.send.Robj(tasks[[1]], slave_id, 1);
        tasks[[1]] <- NULL
      } else {
        mpi.send.Robj(junk, slave_id, 2)
      }
    } else if (tag == 2) {
      response[message[1],] <- message
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
  ret$response <- response

  return(ret)

}
