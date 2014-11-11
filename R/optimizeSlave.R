optimizeSlave <-
function() {
  junk <- 0
  done <- 0

  library(Rmpi)
  
  while (done != 1) {
    # Signal being ready to receive a new task
    mpi.send.Robj(junk,0,1)

    # Receive a task
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    
    if (tag == 1) {
      require(xcms)
      exp_index <- task
      
      if(is.null(xcmsSet_parameters$step)) {     #centWave        
	      print(sapply(xcmsSet_parameters, "[[", exp_index))
		
        xset <- xcmsSet(files=example_sample, method="centWave", 
                  peakwidth=c(xcmsSet_parameters$min_peakwidth[exp_index], xcmsSet_parameters$max_peakwidth[exp_index]),
                  ppm=xcmsSet_parameters$ppm[exp_index], noise=xcmsSet_parameters$noise[exp_index], 
				  snthresh=xcmsSet_parameters$snthresh[exp_index], mzdiff=xcmsSet_parameters$mzdiff[exp_index],
				  prefilter=c(xcmsSet_parameters$prefilter[exp_index], xcmsSet_parameters$value_of_prefilter[exp_index]),
				  mzCenterFun=xcmsSet_parameters$mzCenterFun[exp_index], integrate=xcmsSet_parameters$integrate[exp_index],
				  fitgauss=xcmsSet_parameters$fitgauss[exp_index], verbose.columns=xcmsSet_parameters$verbose.columns[exp_index])
                  
      } else {
      #matchedFilter  
        xset <- xcmsSet(files=example_sample, method="matchedFilter", 
                  fwhm=xcmsSet_parameters$fwhm[exp_index], snthresh=xcmsSet_parameters$snthresh[exp_index],
                  step=xcmsSet_parameters$step[exp_index], steps=xcmsSet_parameters$steps[exp_index],
                  sigma=xcmsSet_parameters$sigma[exp_index], max=xcmsSet_parameters$max[exp_index], 
                  mzdiff=xcmsSet_parameters$mzdiff[exp_index], index=xcmsSet_parameters$index[exp_index])   

      }                   
      
      result <- calcPPS(xset) #, ppm, rt_diff)
      result[1] <- exp_index
   
      rm(xset)
      mpi.send.Robj(result,0,2)
      print("result sent")
    } else if (tag == 2) {
    # Master is saying all tasks are done.  Exit
      done <- 1
    }
    # Else ignore the message or report an error
  }

  # Tell master that this slave is exiting.  Send master an exiting message
  mpi.send.Robj(junk,0,3) 
  
}
