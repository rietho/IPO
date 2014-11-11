optimizeRetGroupSlave <-
function() {
  junk <- 0
  done <- 0
  
  while (done != 1) {
    # Signal being ready to receive a new task
    mpi.send.Robj(junk,0,1)

    # Receive a task
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    print(paste("tag", tag, "task", task))
	
	#print(parameters)
    
    if (tag == 1) {
      require(xcms)

      exp_index <- task

      do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && is.null(parameters$gapInit) && is.null(parameters$gapExtend)
	                && is.null(parameters$plottype) && is.null(parameters$response) 
					&& is.null(parameters$factorDiag) && is.null(parameters$factorGap) && is.null(parameters$localAlignment) 
					&& is.null(parameters$initPenalty)) #&& is.null(parameters$center)) 
					
      retcor_failed = ifelse(do_retcor, 1.1, 1)  
  
      if(do_retcor) {
        try(retcor_failed <- retcor(xset, method="obiwarp", plottype=parameters$plottype[exp_index], distFunc=parameters$distFunc[exp_index],
                             profStep=parameters$profStep[exp_index], center=parameters$center[exp_index], response=parameters$response[exp_index], 
							 gapInit=parameters$gapInit[exp_index], gapExtend=parameters$gapExtend[exp_index],
							 factorDiag=parameters$factorDiag[exp_index], factorGap=parameters$factorGap[exp_index], 
							 localAlignment=parameters$localAlignment[exp_index]))
  	
	
        if(!is.numeric(retcor_failed)) {
          xset <- retcor_failed
          retcor_failed=1
        } 
      }
	  
      minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])      
      try(xset <- group(xset, method="density", bw=parameters$bw[exp_index], mzwid=parameters$mzwid[exp_index], minfrac=minfrac, 
                  minsamp=parameters$minsamp[exp_index], max=parameters$max[exp_index]))	 
 
      result <- calculateRGTV(xset, exp_index, retcor_failed)

      mpi.send.Robj(result,0,2)

    } else if (tag == 2) {
      done <- 1
    }
    # Else ignore the message or report an error
  }

  # Tell master that this slave is exiting.  Send master an exiting message
  mpi.send.Robj(junk,0,3) 
  
}
