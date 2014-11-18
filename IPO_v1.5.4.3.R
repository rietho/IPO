#IPO_V1.3.: good_groups ^ 2 to increase recall 
#IPO_V1.4.: vectorized isotope identification; 
#           no intensity window, between intensity of max carbon and 1
#IPO_V1.5.: in RCSandGSIncreased: also used good_groups ^ 2
#IPO_V1.5.3: no parameter for isotope detection. 
#            c13_peak[,"mz"] has to be within (mzmin + isotope_mass) and (mzmax + isotope_mass)
#            c13_peak[,"rt"] has to be within (rtmin + isotope_mass) and (rtmax + isotope_mass)
#IPO_V1.5.4: if bad_group == 0; bad_group = 1 && good_group += 1
#IPO_V1.5.4.1: changes in calcPPS: 
#                   rt_window <- rt * 0.005
#             		  rt_lower <- part_peaks[,"rt"] - rt_window
#	              	  rt_upper <- part_peaks[,"rt"] + rt_window
#IPO_V1.5.4.2: * added initial parameter check
#              * renamed all factor-variables to params
#              * in optimizeXcmsSet: - also look for mzML-files
#                                    - check if files were found
#			   * bug in optimization for matchedFilter fixed; sigma and mzdiff have to be 
#                    definded later (combineFactors()) when sigma and step as well as steps are already known
#IPO_V1.5.4.3: * LIP calculation in calcPPS fixed 


###methods used for peak picking optimization
###
xcmsSetExperiments <- function(example_sample, params, nSlaves=4) { #ppm=5, rt_diff=0.01, nSlaves=4) {
  library(Rmpi)  
  library(rsm)
    
  junk <- 0
  closed_slaves <- 0
  #nSlaves <- min(mpi.comm.size()-1, nSlaves)  
  
  typ_params <- typeCastFactor(params)
  
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


xcmsSetStatistic <- function(xcms_result, subdir, iterator, score_name="PPS") {

  params <- xcms_result$params
  resp <- xcms_result$response[, "PPS"]
  model <- createModel(xcms_result$design, params$to_optimize, resp)
  xcms_result$model <- model                  
     
  max_settings <- getMaximumExperiment(xcms_result$model)
  tmp <- max_settings[-1]
  tmp[is.na(tmp)] <- 1
  if(!is.null(subdir))
    plotContours(xcms_result$model, tmp, paste(subdir,"/rsm_", iterator, sep=""))

  xcms_result$max_settings <- max_settings

  return(xcms_result)
}

checkXcmsSetParams <- function(params) {
  if(is.null(params$step)) {     #centWave   
    quantitative_parameters <- c("ppm", "min_peakwidth", "max_peakwidth", "snthresh", "mzdiff", "noise", "prefilter", "value_of_prefilter")
    qualitative_parameters <- c("integrate", "fitgauss", "verbose.columns", "mzCenterFun")
    unsupported_parameters <- c("scanrange", "sleep", "ROI.list")
  } else {    
    quantitative_parameters <- c("fwhm", "sigma", "max", "snthresh", "step", "steps", "mzdiff")
    qualitative_parameters <- c("index")
    unsupported_parameters <- c("sleep")  
  } 
  checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
} 


  
checkParams <- function(params, quantitative_parameters, qualitative_parameters, unsupported_parameters) { 
  for(i in 1:length(params)) {
	param <- params[[i]]
	name <- names(params)[i]
	if(name %in% unsupported_parameters) {
	  stop(paste("The parameter", name, "is not supported! Please remove from parameters; stopping!"))
	}
	if(name %in% qualitative_parameters) {
	  if(length(param) == 0) {
	    stop(paste("The parameter", name, "has no value set! Please specify; stopping!"))
	  }
	  if(length(param) > 1) {
	    stop(paste("Optimization of parameter", name, "no supported! Please specify only one value; stopping!"))
	  }
	}
    if(name %in% quantitative_parameters) {
	  if(length(param) == 0) {
	    stop(paste("The parameter", name, "has no value set! Please specify between one and two; stopping!"))
	  }
	  if(length(param) > 2) {
	    stop(paste("Too many values for parameter", name, "! Please specify only one or two; stopping!"))
	  }
	}
  } 
  missing_params <- which(!(c(quantitative_parameters, qualitative_parameters) %in% names(params)))
  if(length(missing_params > 0)) {
   stop(paste("The parameter(s)", paste(c(quantitative_parameters, qualitative_parameters)[missing_params], collapse=", "), "are missing! Please specify; stopping!"))
  }
  
}

optimizeXcmsSet <- function(params=getDefaultStartingXcmsParams(), nSlaves=4, subdir="IPO") { #ppm=5, rt_diff=0.02, nSlaves=4, subdir="IPO") {

  checkXcmsSetParams(params)
   
  example_sample <- list.files(full.names=TRUE, ignore.case=TRUE, recursive=TRUE, pattern="(*.mzX?ML$)") 
  if(length(example_sample)==0)
    stop("no files in directory, stopping!")
	
  centWave <- is.null(params$fwhm)  
  
  history <- list()
  iterator = 1 
  best_range <- 0.25

  if(!file.exists(subdir))
    dir.create(file.path(getwd(), subdir))
    
  while(iterator < 50) {
    cat("\n")
    cat("\n")
    cat("\n")
    cat("starting new DoE with:\n")
    print(params)
        
    xcms_result <- xcmsSetExperiments(example_sample, params, nSlaves) 
#                       ppm, rt_diff, nSlaves)                        
                       
    xcms_result <- xcmsSetStatistic(xcms_result, subdir, iterator)
    history[[iterator]] <- xcms_result     
    params <- xcms_result$params 
    
    if(!resultIncreased(history)) {
      cat("no increase, stopping")
      maxima <- 0
	    max_index <- 1
	    for(i in 1:length(history)) {
	      if(history[[i]]$max_settings[1] > maxima) {
		      maxima <- history[[i]]$max_settings[1]
		      max_index <- i
		    }
	    }
   
        xcms_parameters <- as.list(decodeAll(history[[max_index]]$max_settings[-1], history[[max_index]]$params$to_optimize))      
	    xcms_parameters <- combineParams(xcms_parameters, params$no_optimization)
      
	    if(!is.list(xcms_parameters))
	      xcms_parameters <- as.list(xcms_parameters)
		
	    best_settings <- list()
      best_settings$parameters <- xcms_parameters
      library(xcms)
	  library(Rmpi)
      
      
      if(centWave) {		
        xset <- xcmsSet(files=example_sample, method="centWave", 
                  peakwidth=c(xcms_parameters$min_peakwidth, xcms_parameters$max_peakwidth),
                  ppm=xcms_parameters$ppm, noise=xcms_parameters$noise, 
				          snthresh=xcms_parameters$snthresh, mzdiff=xcms_parameters$mzdiff,
				          prefilter=c(xcms_parameters$prefilter, xcms_parameters$value_of_prefilter),
				          mzCenterFun=xcms_parameters$mzCenterFun, integrate=xcms_parameters$integrate,
				          fitgauss=xcms_parameters$fitgauss, verbose.columns=xcms_parameters$verbose.columns, 
                  nSlaves=nSlaves)
      } else {
        xset <- xcmsSet(files=example_sample, method="matchedFilter", 
                  fwhm=xcms_parameters$fwhm, snthresh=xcms_parameters$snthresh,
                  step=xcms_parameters$step, steps=xcms_parameters$steps,
                  sigma=xcms_parameters$sigma, max=xcms_parameters$max, 
                  mzdiff=xcms_parameters$mzdiff, index=xcms_parameters$index) 
                
                
      }
                
	  best_settings$xset <- xset
      target_value <- calcPPS(xset)#, ppm, rt_diff)
      best_settings$result <- target_value
      history$best_settings <- best_settings
      
      print(best_settings)
      
      return(history)   
    }

         
    for(i in 1:length(params$to_optimize)) {
      parameter_setting <- xcms_result$max_settings[i+1]
      bounds <- params$to_optimize[[i]] 
      fact <- names(params$to_optimize)[i]
	  min_factor <- ifelse(fact=="min_peakwidth", 3, ifelse(fact=="mzdiff", ifelse(centWave,-100000000, 0.001), ifelse(fact=="step",0.0005,1)))

      #if the parameter is NA, we increase the range by 20%, if it was within the inner 25% of the previous range or at the minimum value we decrease the range by 20%
      step_factor <- ifelse(is.na(parameter_setting), 1.2, ifelse((abs(parameter_setting) < best_range), 0.8, ifelse(parameter_setting==-1 & decode(-1, params$to_optimize[[i]]) == min_factor,0.8,1)))
      step <- (diff(bounds) / 2) * step_factor
      
      if(is.na(parameter_setting))
        parameter_setting <- 0
      new_center <- decode(parameter_setting, bounds)
      
      if((new_center-min_factor) > step) {
        new_bounds <- c(new_center - step, new_center + step) 
      } else {
        new_bounds <- c(min_factor, 2*step+min_factor) 
      }      
		  
	  names(new_bounds) <- NULL         
          
      if(names(params$to_optimize)[i] != "mzdiff" & names(params$to_optimize)[i] != "step")
        params$to_optimize[[i]] <- round(new_bounds, 0)
      else 
        params$to_optimize[[i]] <- new_bounds
    } 
    
    if(centWave) {
    #checking peakwidths plausiability
	  if(!is.null(params$to_optimize$min_peakwidth) | !is.null(params$to_optimize$max_peakwidth)) {
	    pw_min <- ifelse(is.null(params$to_optimize$min_peakwidth), params$no_optimization$min_peakwidth, max(params$to_optimize$min_peakwidth))
		  pw_max <- ifelse(is.null(params$to_optimize$max_peakwidth), params$no_optimization$max_peakwidth, min(params$to_optimize$max_peakwidth))
        if(pw_min >= pw_max) {
          additional <- abs(pw_min-pw_max) + 1
          if(!is.null(params$to_optimize$max_peakwidth)) {		  
            params$to_optimize$max_peakwidth <- params$to_optimize$max_peakwidth + additional
          } else {
            params$no_optimization$max_peakwidth <- params$no_optimization$max_peakwidth + additional
		  }
		}
      }
    }
	
	params <- attachparams(params$to_optimize, params$no_optimization)	    
    iterator <- iterator + 1
                 
  }
  return(history)

}


resultIncreased <- function(history) {

  index = length(history)
  if(index < 2)
    return(TRUE)
   
  if(history[[index-1]]$max_settings[1] >= history[[index]]$max_settings[1])
    return(FALSE)
    
  return(TRUE)

}

getDefaultStartingXcmsParams <- function(method="centWave") {
  if(method=="centWave")
    return(list(min_peakwidth=c(10,30), max_peakwidth=c(35,65), ppm=c(15,35),
              mzdiff=c(-0.001, 0.01), snthresh=10, noise=0, prefilter=3, 
			        value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
			        fitgauss=FALSE, verbose.columns=FALSE))
			  
  if(method=="matchedFilter")
    return(list(fwhm=c(25,35), snthresh=c(3,17), step=c(0.05, 0.15), steps=c(1,3), 
           sigma=0, max=5, mzdiff=0, index=FALSE)) 
  
}

startSlaves <- function(nSlaves) {
  mpi.spawn.Rslaves(nslaves=nSlaves)
                                                                                
  .Last <- function() {
    if (is.loaded("mpi_initialize")) {
      if (mpi.comm.size(1) > 0) {
	    print("Please use mpi.close.Rslaves() to close slaves.")
	    mpi.close.Rslaves()
	  }
      print("Please use mpi.quit() to quit R")
      .Call("mpi_finalize")
    }
  }
}

sendXcmsSetSlaveFunctions <- function(example_sample, xcmsSet_parameters) {#, ppm, rt_diff) {
  mpi.bcast.Robj2slave(toMatrix) 
  mpi.bcast.Robj2slave(calcPPS) 
  mpi.bcast.cmd(slave <- mpi.comm.rank())
  mpi.bcast.Robj2slave(example_sample)
  mpi.bcast.Robj2slave(xcmsSet_parameters)
  mpi.bcast.Robj2slave(optimizeSlave)
  
  mpi.bcast.cmd("library(xcms)")
  mpi.bcast.cmd(optimizeSlave())
}

optimizeSlave <- function() {
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

calcPPS <- function(xset) { 
   
  iso_list <- list()
  ret <- array(0, dim=c(1,5))  
  if(is.null(xset)) {
    return(ret)
  }
  
  peak_source <- toMatrix(xset@peaks[,c("mz", "rt", "sample", "into", "mzmin", "mzmax", "rtmin", "rtmax")])
  if(nrow(peak_source) == 0) {
	  return(ret)
  }
  
  for(i in 1:ncol(peak_source)) {
    peak_source <- toMatrix(peak_source[!is.na(peak_source[,i]),])
  }
  
  ret[2] <- nrow(peak_source)

  peak_source <- cbind(1:nrow(peak_source), peak_source)
  colnames(peak_source)[1] <- "id"  
  
  carbon = 12.0
  hydrogen	= 1.0078250170
  CH3 = carbon + 3 * hydrogen
  CH2 = carbon + 2 * hydrogen
  isotope_mass = 1.0033548
  samples <- max(peak_source[,"sample"])

  #start_sample
  for(sample in 1:samples) { 
    #only taking peaks from current sample   
	  speaks <- toMatrix(peak_source[peak_source[,"sample"]==sample,])	
	  found_isotope <- FALSE
    split <- 250
    	
	  if(nrow(speaks)>1) {  		      
	    #speaks <- speaks[,-c("sample")]
	    speaks <- speaks[order(speaks[,"mz"]),]
		      
	    while(!is.null(nrow(speaks)) & length(speaks) > 3) {
	      part_peaks <- NULL
	      #splitting the data into smaller pieces to improve speed    
	      if(nrow(speaks) < split) {
	        part_peaks <- speaks
	      } else {          
            upper_bound <- speaks[split,"mzmax"] + isotope_mass# + (speaks[split,"mz"] + isotope_mass) * ppm / 1000000          
	        end_point <- sum(speaks[,"mz"] < upper_bound)
	        part_peaks <- toMatrix(speaks[1:end_point,])
	      }		

		    rt <- part_peaks[,"rt"]
		    rt_window <- rt * 0.005
		    rt_lower <- part_peaks[,"rt"] - rt_window
		    rt_upper <- part_peaks[,"rt"] + rt_window
		    rt_matrix <-  t(matrix(rep(rt, nrow(part_peaks)), ncol=nrow(part_peaks)))
		    rt_matrix_bool <- rt_matrix >= rt_lower & rt_matrix <= rt_upper

		    mz <- part_peaks[,"mz"]
		    mz_lower <- part_peaks[,"mzmin"] + isotope_mass #isotope_masses - mz_window
		    mz_upper <- part_peaks[,"mzmax"] + isotope_mass #isotope_masses + mz_window
		    mz_matrix <-  t(matrix(rep(mz, nrow(part_peaks)), ncol=nrow(part_peaks)))
		    mz_matrix_bool <- mz_matrix >= mz_lower & mz_matrix <= mz_upper

		    rt_mz_matrix_bool <- rt_matrix_bool & mz_matrix_bool
		  
		    rt_mz_peak_ids <- which(rowSums(rt_mz_matrix_bool)>0)
		    calculations <- min(split, nrow(speaks))
		    rt_mz_peak_ids <- rt_mz_peak_ids[rt_mz_peak_ids < calculations]
      
        #if(length(rt_mz_peak_ids)>0) {
		      for(i in rt_mz_peak_ids) {
			      current <- part_peaks[i,]
			      rt_mz_peaks <- toMatrix(part_peaks[rt_mz_matrix_bool[i,],])
			      rt_difference <- abs(current["rt"] - rt_mz_peaks[, "rt"]) / current["rt"]
			      rt_mz_peaks <- cbind(rt_mz_peaks, rt_difference)
            #test intensity_window
            maximum_carbon <- floor((current["mz"]-2*CH3)/CH2) + 2
            carbon_probabilty <- c(1,maximum_carbon)*0.01108
            iso_intensity <- current["into"] * carbon_probabilty

            int_bools <- rt_mz_peaks[,"into"] >= iso_intensity[1] & rt_mz_peaks[,"into"] <= iso_intensity[2]
            if(sum(int_bools) > 0) {
              int_peaks <- toMatrix(rt_mz_peaks[int_bools,])
              iso_id <- int_peaks[which.min(int_peaks[,"rt_difference"]), "id"]
              iso_list[[length(iso_list)+1]] <- c(current["id"], iso_id)  
              found_isotope <- TRUE              
            }
		      }
        #}
		    speaks <- speaks[-(1:calculations),]		    
	      
      }#end_while_sample_peaks      

      sample_isos_peaks <- peak_source
      sample_non_isos_peaks <- peak_source
      
      if(found_isotope) {
        sample_isos_peaks <- toMatrix(peak_source[unique(unlist(iso_list)),])
        sample_non_isos_peaks <- toMatrix(peak_source[-unique(unlist(iso_list)),])
      } 

      speaks <- toMatrix(sample_non_isos_peaks[sample_non_isos_peaks[,"sample"]==sample,])
      sample_isos_peaks <- toMatrix(sample_isos_peaks[sample_isos_peaks[,"sample"]==sample,])
      int_cutoff = 0
      iso_int <- speaks[,"into"]

      tmp <- iso_int[order(iso_int)]      
      int_cutoff <- mean(tmp[1:round((length(tmp)/33),0)])

      masses <- speaks[, "mz"]
      maximum_carbon <- floor((masses-2*CH3)/CH2) + 2
      carbon_probabilty <- maximum_carbon*0.01108

      iso_int <- iso_int * carbon_probabilty
  
      not_loq_peaks <- sum(iso_int>int_cutoff)
      ret[3] <- ret[3] + not_loq_peaks
      ret[4] <- length(unique(unlist(iso_list)))
	  if(ret[3] == 0) {
	    ret[5] <- (ret[4]+1)^1.5/(ret[3]+1)  
	  } else {	  
        ret[5] <- ret[4]^1.5/ret[3]  
	  }
	  }    
  }#end_for_sample    
  
  return(ret)

}


xcmsSetsettingsAsString <- function(parameters) {  
  
  return(sprintf("xcmsSet(method='centWave', ppm=%s, peakwidth=c(%s,%s), snthresh=%s, prefilter=c(%s,%s), mzCenterFun='%s', integrate=%s, mzdiff=%s, fitgauss=%s, noise=%s, verbose.columns=%s)", 
          prettyNum(parameters$ppm), prettyNum(parameters$min_peakwidth), 
          prettyNum(parameters$max_peakwidth), prettyNum(parameters$snthresh), 
          prettyNum(parameters$prefilter), prettyNum(parameters$value_of_prefilter), 
          parameters$mzCenterFun, prettyNum(parameters$integrate), 
          prettyNum(parameters$mzdiff), parameters$fitgauss, 
          prettyNum(parameters$noise), parameters$verbose.columns))
}

###methods used for retention time correction and grouping optimization
###
checkRetGroupSetParams <- function(params) {

  quantitative_parameters <- c("profStep", "gapInit", "gapExtend", "response", "factorDiag", "factorGap", "minfrac", "minsamp", "bw", "mzwid", "max")
  qualitative_parameters <- c("distFunc", "plottype", "localAlignment", "center")
  unsupported_parameters <- c("col", "ty", "initPenalty", "sleep")

  checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
} 

retGroupCalcExperiments <- function(params, xset, nSlaves=4) {
							   
  library(Rmpi)  
  library(rsm)
    
  junk <- 0
  closed_slaves <- 0
  #nSlaves <- min(mpi.comm.size()-1, nSlaves)  
  
  typ_params <- typeCastFactor(params)
  
  if(length(typ_params$to_optimize) > 2) {
    design <- getBbdParameter(typ_params$to_optimize) 
  } else {
    design <- getCcdParameter(typ_params$to_optimize) 
  }	
 
  parameters <- decode.data(design)	
  tasks <- as.list(1:nrow(design))    
  startSlaves(nSlaves)
  
  parameters <- combineParams(parameters, typ_params$no_optimization)
  
  sendRetGroupSlaveFunctions(parameters, xset) 
  
  response <- list()  
  finished <- 0
  while(closed_slaves < nSlaves) {
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

getDefaultRetGroupStartingParams <- function(distfunc="cor_opt", high_resolution=TRUE) {

  ret <- NULL
  if(!is.null(distfunc)) {
    if(distfunc=="cor")
      ret <- (list(distFunc="cor", gapInit=c(0.0, 0.4), gapExtend=c(2.1, 2.7)))
	if(distfunc=="cor_opt")
	  ret <- (list(distFunc="cor_opt", gapInit=c(0.0, 0.4), gapExtend=c(2.1, 2.7)))
    if(distfunc=="cov")
	  ret <- (list(distFunc="cov", gapInit=c(0.0, 0.4), gapExtend=c(11.4, 12.0)))
	if(distfunc=="prd")
	  ret <- (list(distFunc="prd", gapInit=c(0.0, 0.4), gapExtend=c(7.5, 8.1)))
	if(distfunc=="euc")
	  ret <- (list(distFunc="euc", gapInit=c(0.7, 1.1), gapExtend=c(1.5, 2.1)))

    ret$profStep <- c(0.7, 1)
    ret$plottype <- "none"
    ret$response <- 1
    ret$factorDiag <- 2
    ret$factorGap <- 1
    ret$localAlignment <- 0
    #ret$initPenalty <- 0
	  
  } else {
	ret <- list()
  }

	#grouping parameter
  ret$bw <- c(22,38)
  ret$minfrac <- c(0.3, 0.7)
  ret$mzwid <- c(0.15, 0.35)
  if(high_resolution)
    ret$mzwid <- c(0.015, 0.035)
  ret$minsamp <- 1
  ret$max <- 50
  #ret$sleep <- 0
	  
  return(ret)

}

calculateRGTV <- function(xset, exp_index=1, retcor_penalty=1) {
  features <- nrow(xset@groups)

  good_groups <- sum(unlist(lapply(X=xset@groupidx, FUN=function(x, xset) {ifelse(length(unique(xset@peaks[x,"sample"]))==length(xset@filepaths) & 
                                   length(xset@peaks[x,"sample"])==length(xset@filepaths),1,0)}, xset)))
  bad_groups <- nrow(xset@groups) - good_groups

  relative_rt_diff <- c()
  
  if(nrow(xset@groups) > 0) {
    for(i in 1:nrow(xset@groups)) {
      feature_rtmed <- xset@groups[i, "rtmed"]
	    relative_rt_diff <- c(relative_rt_diff, mean(abs(feature_rtmed - xset@peaks[xset@groupidx[[i]], "rt"])/feature_rtmed))
    }
  } else {
    relative_rt_diff <- 1
  }
  
  ARTS <- (mean(relative_rt_diff)) * retcor_penalty
  
  ret <- list(exp_index=exp_index, good_groups=good_groups, bad_groups=bad_groups, mean_rel_rt_diff=ARTS)
  
  ret$retcor_done = retcor_penalty        
  
  return(ret)  
}

retGroupExperimentStatistic<- function(retcor_result, subdir, iterator, xset) {

  params <- retcor_result$params
  resp <- getNormalizedResponse(retcor_result$response)
  model <- createModel(retcor_result$design, params$to_optimize, resp)
  
  retcor_result$model <- model                  
  max_settings <- getMaximumExperiment(retcor_result$model)
  tmp <- max_settings[-1]
  tmp[is.na(tmp)] <- 1

  if(!is.null(subdir))
    plotContours(retcor_result$model, tmp, paste(subdir, "/retgroup_rsm_", iterator, sep=""))  
    
  parameters <- as.list(decodeAll(max_settings[-1], params$to_optimize)) 
  parameters <- combineParams(parameters, params$no_optimization)
  xset_tmp <- xset
  
  exp_index <- 1				   
  do_retcor <- !(is.null(parameters$distFunc) && is.null(parameters$profStep) && is.null(parameters$gapInit) && is.null(parameters$gapExtend)
	             && is.null(parameters$plottype) && is.null(parameters$col) && is.null(parameters$ty) && is.null(parameters$response) 
				 && is.null(parameters$factorDiag) && is.null(parameters$factorGap) && is.null(parameters$localAlignment) 
				 && is.null(parameters$initPenalty)) #&& is.null(parameters$center)) 
					
      retcor_failed = ifelse(do_retcor, 1.1, 1)  
  
  if(do_retcor) {
    try(retcor_failed <- retcor(xset_tmp, method="obiwarp", plottype=parameters$plottype[exp_index], distFunc=parameters$distFunc[exp_index],
                             profStep=parameters$profStep[exp_index], center=parameters$center[exp_index],
							 response=parameters$response[exp_index], gapInit=parameters$gapInit[exp_index], 
				 	         gapExtend=parameters$gapExtend[exp_index], factorDiag=parameters$factorDiag[exp_index], 
					         factorGap=parameters$factorGap[exp_index], localAlignment=parameters$localAlignment[exp_index]))
  	
	
    if(!is.numeric(retcor_failed)) {
      xset_tmp <- retcor_failed
      retcor_failed=1
    } 
  }
	  
  minfrac <- ifelse(is.null(parameters$minfrac), 1, parameters$minfrac[exp_index])      
  try(xset_tmp <- group(xset_tmp, method="density", bw=parameters$bw[exp_index], mzwid=parameters$mzwid[exp_index], minfrac=minfrac, 
                  minsamp=parameters$minsamp[exp_index], max=parameters$max[exp_index]))	 
 
  tv <- calculateRGTV(xset_tmp, exp_index, retcor_failed)

  retcor_result$max_settings <- max_settings
  retcor_result$target_value <- tv   
  retcor_result$xset <- xset_tmp
  retcor_result$best_settings <- parameters

  return(retcor_result)

}

getNormalizedResponse <- function(response) {

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

optimizeRetGroup <- function(xset, params=getDefaultRetGroupStartingParams(), nSlaves=4, subdir="IPO") {
                                                 
  library(xcms)
  iterator = 1 
  history <- list()  
  best_range <- 0.25

  if(!is.null(subdir) & !file.exists(subdir))
    dir.create(file.path(getwd(), subdir))
	
  if(is.null(params$center))
    params$center <- getDefaultRetCorCenterSample(xset) 

  checkRetGroupSetParams(params)	
    
  while(iterator < 50) {
    cat("\n")
    cat("\n")
    cat("\n")
    cat("starting new DoE with:\n")
    print(params)
        
    retcor_result <- retGroupCalcExperiments(params, xset, nSlaves)  
                       
    retcor_result <- retGroupExperimentStatistic(retcor_result, subdir, iterator, xset)
    
    history[[iterator]] <- retcor_result 
         
    params <- retcor_result$params$to_optimize 
	#all settings must be within range, not like xcmsSet-Optimization
    #found_good_settings <- foundGoodSettings(rep(TRUE, length(retcor_result$parameter_significance)), retcor_result$eigenvalues, retcor_result$max_settings[-1])
    
    if(!RCSandGSIncreased(history)) {
      print("no increase stopping")
      
      history$best_settings <- history[[(length(history)-1)]]$best_settings	 
	  
      return(history)       

    }
             
    for(i in 1:length(params)) {

      parameter_setting <- retcor_result$max_settings[i+1]
      bounds <- params[[i]] 
		
	    min_bound <- ifelse(names(params)[i] == "profStep",0.3, ifelse(names(params)[i] == "mzwid",0.0001,ifelse(names(params)[i] == "bw",0.25,0)))		
      
      step_factor <- ifelse(is.na(parameter_setting), 1.2, ifelse((abs(parameter_setting) < best_range), 0.8, 1))
      step <- (diff(bounds) / 2) * step_factor
      
      if(is.na(parameter_setting)) {      
		    new_center <- decode(0, bounds)
		    step <- (diff(bounds) / 2) * 1.2
		    if(new_center > step) {
		      new_bounds <- c(new_center - step, new_center + step) 
		    } else {
		      new_bounds <- c(min_bound, 2*step + min_bound)      
		    }

		    names(new_bounds) <- NULL
            params[[i]] <- new_bounds
			
		  } else {		  
		    if(bounds[1] > min_bound | parameter_setting > -1) {
          #center around optimum
          new_center <- decode(parameter_setting, bounds)
          step <- diff(bounds) / 2
	      if(abs(parameter_setting) < best_range) { #if parameter are within range, decrease parameter-range
	        step <- step * 0.8
        }

			  if(new_center > step) {
          new_bounds <- c(new_center - step, new_center + step) 
			  } else {
		      new_bounds <- c(min_bound, 2*step + min_bound)  
        }			  

		    if(parameter_setting == -1 || parameter_setting == 1) {
		      step <- step / 5
		      new_bounds[1] <- max(min_bound, new_bounds[1] - step)
		      new_bounds[2] <- new_bounds[2] + step
		    } 
			#print(paste("new_bounds", paste(new_bounds, collapse=" ")))
            names(new_bounds) <- NULL
			
            params[[i]] <- new_bounds
          } else {
		        step <- diff(bounds) / 2
            params[[i]] <- c(min_bound, 2*step + min_bound)
		      }		  
      }
      
      if(names(params)[i] == "profStep" | names(params)[i] == "minfrac") {
        if(params[[i]][2] > 1) {  # 1 is max value for profStep
          params[[i]] <- round(c(1-(diff(params[[i]])*0.8), 1),2)
          print(paste("profStep or minfrac greater 1, decreasing to", paste(params[[i]], collapse = " and ")))  
        }      
      }
		
	}
	
	params <- attachparams(params, retcor_result$params$no_optimization)
              
    iterator <- iterator + 1
                 
  }
  return(history)

}

sendRetGroupSlaveFunctions <- function(parameters, xset) {
  mpi.bcast.Robj2slave(calculateRGTV)
  mpi.bcast.cmd(slave <- mpi.comm.rank())
  mpi.bcast.Robj2slave(parameters)
  mpi.bcast.Robj2slave(xset)  
  mpi.bcast.cmd("library(xcms)")
  mpi.bcast.Robj2slave(optimizeRetGroupSlave)  
  mpi.bcast.cmd(optimizeRetGroupSlave())
}

optimizeRetGroupSlave <- function() {
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

RCSandGSIncreased <- function(history) {

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

###general methods used for optimization of peak picking and retention time correction and grouping
###
combineParams <- function(params_1, params_2) {
  len <- max(unlist(sapply(params_1, length)))
  #num_params <- length(params_1)
  
  p_names <- c(names(params_1), names(params_2))
  matchedFilter <- "fwhm" %in% p_names
  
  for(i in 1:length(params_2)) {
    new_index <- length(params_1) + 1
	fact <- params_2[[i]]
	params_1[[new_index]] <- fact
	if(matchedFilter) {
	  if(p_names[new_index] == "sigma" && fact == 0) { #update values for sigma if zero
	    if("fwhm" %in% names(params_1)) {
	      params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
	    } else {
	      params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
	    }	
	  } else {
	    if(p_names[new_index] == "mzdiff" && fact == 0) { #update values for mzdiff if zero		
		  if("step" %in% names(params_1)) {
	        if("steps"  %in% names(params_1)) {
	          params_1[[new_index]][1:len] <- 0.8-params_1$step*params_1$steps
	        } else {
	          params_1[[new_index]][1:len] <- 0.8-params_1$step*params_2$steps
	        }	
	      } else {
	        if("steps"  %in% names(params_1)) {
	          params_1[[new_index]][1:len] <- 0.8-params_2$step*params_1$steps
	        } else {
	          params_1[[new_index]][1:len] <- 0.8-params_2$step*params_2$steps
	        }	
	      }		  
		} else {          
          params_1[[new_index]][1:len] <- fact
		}
	  }
	} else {
      params_1[[new_index]][1:len] <- fact
	}
  } 

  names(params_1) <- p_names   
  return(params_1)

}

attachparams <- function(params_1, params_2) {
  params <- params_1
  for(factor in params_2)
	params[[length(params)+1]] <- factor
	  
  names(params) <- c(names(params_1), names(params_2))
  return(params)
}

#split parameters in those which should be optimized and which should not
typeCastFactor <- function(params) {
  ret_1 <- list()
  ret_2 <- list()  
  ret <- list()
  for(i in  1:length(params)) {
    factor <- params[[i]]
    if(length(factor) == 2) {
	  ret_1[[(length(ret_1)+1)]] <- factor
	  names(ret_1)[length(ret_1)] <- names(params)[i]
	} else {	
	  ret_2[[(length(ret_2)+1)]] <- factor
	  names(ret_2)[length(ret_2)] <- names(params)[i]
	}
  }	
  ret$to_optimize <- ret_1
  ret$no_optimization <- ret_2
  
  return(ret)
}

getMaximumExperiment <- function(model) {  
  dimensions <- length(model$coding)   
  slices <- getSlices(dimensions-2)
  mat <- getResponses(slices, model)  
  return(mat[which.max(mat[,1]),])
}

encode <- function(value, bounds) {
  x <- (value - min(bounds)) / (max(bounds) - min(bounds))
  
  return(x*2-1)
}

decode <- function(value, bounds) {
  if(is.na(value))
    value <- 1
  x <- (value+1)/2
  x <- (x*(max(bounds)-min(bounds))) + min(bounds)
  
  return(x)
}

decodeAll <- function(values, params) {

  ret <- rep(0, length(params))
  for(i in 1:length(params))
    ret[i] <- decode(values[i], params[[i]])
  
  names(ret) <- names(params)
  
  return(ret)
}

getResponses <- function(slices, model) {

   if(is.null(slices)) {
     slices= matrix(nrow=1, ncol=0)
   }
   add_col <- array(0, dim=c(nrow(slices), 2))
   slices <- cbind(add_col, slices)
   
   colnames(slices) <- paste("x", 1:ncol(slices), sep="")
   values <- apply(X=slices, MARGIN=1, FUN=getMaxima, model) 
    
   ret <- cbind(t(values), slices[,-c(1,2)]) 
   colnames(ret) <- c("response", paste("x", 1:ncol(slices), sep=""))
    
   return(ret)
}

getMaxima <- function(slice, model) {
  	
  values <- contour(model, ~ x1 * x2, at=slice, plot.it=FALSE, decode=FALSE) 
    
  z_values <- values[[1]]$z
  max_z <- max(z_values)
  maxima <- which(z_values==max_z, arr.ind=T)
  ret <- c(max_z, 0,0)
  if(length(maxima)==2){
    ret[2] <- values[[1]]$x[maxima[1]]
    ret[3] <- values[[1]]$y[maxima[2]]     
  } else {
    for(i in 1:2) {
      col_max <- maxima[,i]
	    if(sum(col_max==1) > 0 && sum(col_max==length(values[[1]][[i]])) > 0) 
	      ret[(i+1)] <- NA
      else
   	    ret[(i+1)] <- values[[1]][[i]][max(col_max)]
	  }
  }
  
  return(ret)
}

getSlices <- function(dimensions, slice=NULL) {
  ret <- c()
  if(length(slice) == dimensions) {
    return(slice)
  } else {
    values <- seq(-1,1,0.2)
    for(val in values)
      ret <- rbind(ret, getSlices(dimensions, c(slice, val)))
  }
  
  return(ret)
}

plotContours <- function(model, maximum_slice, plot_name) {

  plots <- c()
  for(i in 1:(length(maximum_slice)-1)) {
    for(j in (i+1):length(maximum_slice)) {
      plots <- c(plots, as.formula(paste("~ x", i, "* x", j, sep="")))
    } 
  }
    
  plot_rows <- round(sqrt(length(plots)))
  plot_cols <- if(plot_rows==1){length(plots)}else{ceiling(sqrt(length(plots)))}

  if(!is.null(plot_name)) {
    plot_name = paste(plot_name, ".jpg", sep="")
    jpeg(plot_name, width=4*plot_cols, height=2*plot_rows+2, unit="in", res=c(200,200))
  } else 
    dev.new(width=4*plot_cols, height=2*plot_rows+2)
  par(mfrow=c(plot_rows, plot_cols), oma=c(3,0,2,0))  
  contours <- contour(model, plots, image=TRUE, at=maximum_slice)
  
  if(!is.null(plot_name)) {
    dev.off()
  }

}

createModel <- function(design, params, resp) {
  design$resp <- resp
  formula <- as.formula(paste("resp ~ SO(", paste("x", 1:length(params), sep="", collapse=","), ")", sep="")) 
  return(rsm(formula, data=design)) 
}

getCcdParameter <- function(params) {
 
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

getBbdParameter <- function(params) {
 
  require(rsm)
  lower_bounds <- unlist(lapply(X=params, FUN=function(x) x[1]))
  higher_bounds <- unlist(lapply(X=params, FUN=function(x) x[2]))
  
  steps <- (higher_bounds - lower_bounds)/2

  x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", (lower_bounds + steps), ")/", steps, sep="")
  formulae <- list()
  for(i in 1:length(x))
    formulae[[i]] <- as.formula(x[i])  
  
  design <- bbd(length(params), n0 = 1, randomize = FALSE, coding = formulae)

  return(design)
  
}

toMatrix <- function(data) {

  if(!is.matrix(data)) {
    tmp <- names(data)
    data <- matrix(data, nrow=1)
    colnames(data) <- tmp  
  } 
  
  return(data)  
  
}

getDefaultRetCorCenterSample <- function(xset) {
  ret <- NULL
  for(i in 1:length(xset@filepaths)) {
    ret <- c(ret, sum(xset@peaks[,"sample"] == i))
  }
  return(which.max(ret))
}

writeRSkript <- function(peakPickingSettings, retCorGroupSettings, nSlaves) {
  cat("library(xcms)\n")
  cat("library(Rmpi)\n")

  if(is.null(peakPickingSettings$step)) {     #centWave     		
    cat(paste("xset <- xcmsSet(method=\"centWave\", peakwidth=c(", 
              peakPickingSettings$min_peakwidth, ", ", peakPickingSettings$max_peakwidth,
              "), ppm=", peakPickingSettings$ppm, ", noise=", peakPickingSettings$noise, 
			  ", snthresh=", peakPickingSettings$snthresh, ", mzdiff=", peakPickingSettings$mzdiff,
			  ", prefilter=c(", peakPickingSettings$prefilter, ", ", peakPickingSettings$value_of_prefilter,
			  "), mzCenterFun=\"", peakPickingSettings$mzCenterFun, "\", integrate=", peakPickingSettings$integrate,
			  ", fitgauss=", peakPickingSettings$fitgauss, ", verbose.columns=", peakPickingSettings$verbose.columns,
			  ", nSlaves=", nSlaves, ")\n", sep=""))
                  
  } else { #matchedFilter  
    cat(paste("xset <- xcmsSet(method=\"matchedFilter\", fwhm=", 
	         peakPickingSettings$fwhm, ", snthresh=",peakPickingSettings$snthresh,
             ", step=", peakPickingSettings$step, ", steps=", peakPickingSettings$steps,
             ", sigma=", peakPickingSettings$sigma, ", max=", peakPickingSettings$max, 
             ", mzdiff=", peakPickingSettings$mzdiff, ", index=", peakPickingSettings$index,
	         ")\n", sep=""))   
  }
	  
	  
  cat(paste("xset <- retcor(xset, method=\"obiwarp\", plottype=\"", retCorGroupSettings$plottype, 
            "\", distFunc=\"", retCorGroupSettings$distFunc, "\", profStep=", retCorGroupSettings$profStep, 
			", center=", retCorGroupSettings$center, ", response=", retCorGroupSettings$response, 
		    ", gapInit=", retCorGroupSettings$gapInit, ", gapExtend=", retCorGroupSettings$gapExtend,
		    ", factorDiag=", retCorGroupSettings$factorDiag, ", factorGap=", retCorGroupSettings$factorGap, 
			", localAlignment=", retCorGroupSettings$localAlignment, ")\n", sep=""))
  	   
  cat(paste("xset <- group(xset, method=\"density\", bw=", retCorGroupSettings$bw, 
            ", mzwid=", retCorGroupSettings$mzwid, ", minfrac=", retCorGroupSettings$minfrac, 
            ", minsamp=", retCorGroupSettings$minsamp, ", max=", retCorGroupSettings$max,
	        ")\n", sep=""))	 
	  
  cat(paste("xset <- fillPeaks(xset, nSlaves=", nSlaves, ")\n", sep=""))	
	
}


#setwd("/...")

#peakpicking_parameter <- getDefaultStartingXcmsParams()
#result_peakpicking <- optimizeXcmsSet(params=peakpicking_parameter, nSlaves=4)

#retcor_parameter <- getDefaultRetGroupStartingParams()
#result_retcor <- optimizeRetGroup(result_peakpicking$best_settings$xset, params=retcor_parameter, nSlaves=4)


