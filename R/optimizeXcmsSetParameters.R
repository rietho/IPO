calcPPS <-
  function(xset, isotopeIdentification=c("IPO", "CAMERA"), ...) {
    
    isotopeIdentification <- match.arg(isotopeIdentification)
    
    ret <- vector(mode="numeric", 5) #array(0, dim=c(1,5)) 
    names(ret) <- c("ExpId", "#peaks", "#NonRP", "#RP", "PPS")
    if(is.null(xset)) {
      return(ret)
    } 
    
    if(nrow(xset@peaks) == 0) {
      return(ret)
    }
    
    peak_source <- xset@peaks[,c("mz", "rt", "sample", "into", "mzmin", "mzmax", "rtmin", "rtmax"),drop=FALSE]
    ret[2] <- nrow(peak_source)
    
    if(isotopeIdentification == "IPO")
      iso_mat <- findIsotopes.IPO(xset, ...)  
    else
      iso_mat <- findIsotopes.CAMERA(xset, ...)
    
    samples <- unique(peak_source[,"sample"])
    isotope_abundance = 0.01108    
    
    #calculating low intensity peaks
    for(sample in samples) {
      non_isos_peaks <- peak_source
      
      if(nrow(iso_mat) > 0) {
        non_isos_peaks <- toMatrix(peak_source[-unique(c(iso_mat)),])
      } 
      
      speaks <- toMatrix(non_isos_peaks[non_isos_peaks[,"sample"]==sample,])
      iso_int <- speaks[,"into"]
      
      tmp <- iso_int[order(iso_int)]      
      int_cutoff <- mean(tmp[1:round((length(tmp)/33),0)])
      
      masses <- speaks[, "mz"]
      maximum_carbon <- calcMaximumCarbon(masses)#floor((masses-2*CH3)/CH2) + 2
      carbon_probabilty <- maximum_carbon*isotope_abundance
      
      iso_int <- iso_int * carbon_probabilty
      
      not_loq_peaks <- sum(iso_int>int_cutoff)
      ret[3] <- ret[3] + not_loq_peaks
      
    }#end_for_sample    
    
    ret[4] <- length(unique(c(iso_mat)))
    if(ret[3] == 0) {
      ret[5] <- (ret[4]+1)^2/(ret[3]+1)  
    } else {    
      ret[5] <- ret[4]^2/ret[3]  
    }
    
    return(ret)
    
  }




findIsotopes.IPO <- 
  function(xset, checkBorderIntensity=FALSE) {
    
    iso_mat <- matrix(0, nrow=0, ncol=2)
    colnames(iso_mat) <- c("12C", "13C")
    peak_source <- toMatrix(xset@peaks[,c("mz", "rt", "sample", "into", "maxo", "mzmin", "mzmax", "rtmin", "rtmax")])
    
    for(i in 1:ncol(peak_source)) {
      peak_source <- toMatrix(peak_source[!is.na(peak_source[,i]),])
    }
    
    peak_source <- cbind(1:nrow(peak_source), peak_source)
    colnames(peak_source)[1] <- "id"  
    
    #carbon = 12.0
    #hydrogen	= 1.0078250170
    #CH3 = carbon + 3 * hydrogen
    #CH2 = carbon + 2 * hydrogen
    isotope_mass = 1.0033548
    isotope_abundance = 0.01108
    
    samples <- max(peak_source[,"sample"])
    
    #start_sample
    for(sample in 1:samples) { 
      #only looking into peaks from current sample   
      speaks <- toMatrix(peak_source[peak_source[,"sample"]==sample,])	
      split <- 250
      rawdata <- loadRaw(xcmsSource(xset@filepaths[sample]))
      
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
          
          for(i in rt_mz_peak_ids) {
            current <- part_peaks[i, ,drop=FALSE]
            rt_mz_peaks <- part_peaks[rt_mz_matrix_bool[i,],,drop=FALSE]
            rt_difference <- abs(current[,"rt"] - rt_mz_peaks[, "rt"]) / current[,"rt"]
            rt_mz_peaks <- cbind(rt_mz_peaks, rt_difference)
            #test intensity_window
            maximum_carbon <- calcMaximumCarbon(current[,"mz"]) #floor((current["mz"]-2*CH3)/CH2) + 2
            carbon_probabilty <- c(1,maximum_carbon)*isotope_abundance
            iso_intensity <- current[,"into"] * carbon_probabilty
            
            int_bools <- rt_mz_peaks[,"into"] >= iso_intensity[1] & rt_mz_peaks[,"into"] <= iso_intensity[2]
            
            if(sum(int_bools) > 0) {
              int_peaks <- rt_mz_peaks[int_bools,,drop=FALSE]
              if(checkBorderIntensity) {
                boundary_bool <- checkIntensitiesAtRtBoundaries(rawdata, rbind(current,int_peaks[,-ncol(int_peaks)]))
              } else {
                boundary_bool <- rep(TRUE, (nrow(int_peaks)+1))
              }              
              if(boundary_bool[1] & sum(boundary_bool[-1])>0) {                 
                iso_peaks <- int_peaks[boundary_bool[-1],,drop=FALSE]
                iso_id <- iso_peaks[which.min(iso_peaks[,"rt_difference"]), "id"]
                #iso_list[[length(iso_list)+1]] <- c(current[,"id"], iso_id)            
                iso_mat <- rbind(iso_mat, c(current[,"id"], iso_id))
                
              }
            }
          }
          speaks <- speaks[-(1:calculations),]		    
          
        }#end_while_sample_peaks 
      }
    }
    return(iso_mat)
  }

#checking intensities at rtmin and rtmax. peaks[,"maxo"] must be at least double as high
#does not work for retention time corrected data 
checkIntensitiesAtRtBoundaries <- function(rawdata, peaks, minBoundaryToMaxo=1/3, ppmstep=15) {
  ret <- rep(TRUE, nrow(peaks))
  for(i in 1:nrow(peaks)) {
    peak <- peaks[i,]
    for(boundary in c("rtmin", "rtmax")) {
      rtIndex <- which(rawdata$rt==peak[boundary])
      if(length(rtIndex)>0) {
        if(rtIndex==length(rawdata$scanindex)) {
          rtIndices <- c(rawdata$scanindex[rtIndex], length(rawdata$mz))
        } else {
          rtIndices <- rawdata$scanindex[c(rtIndex, rtIndex+1)]
        }
        
        #only relevant mz and intensity values regarding retention time
        mz <- rawdata$mz[(rtIndices[1]+1):rtIndices[2]]	
        intensities <- rawdata$intensity[(rtIndices[1]+1):rtIndices[2]]
        
        ppm <- peak[c("mzmin", "mzmax")]*ppmstep/1000000
        mzIntensities <- c(0,intensities[mz>=peak["mzmin"]-ppm[1] & mz<=peak["mzmax"]+ppm[2]])
        maxBoundaryIntensity <- max(mzIntensities)
        ret[i] <- ret[i] & maxBoundaryIntensity<peak["maxo"]*minBoundaryToMaxo
      }
    }
  }
  
  return(ret)
  
}

findIsotopes.CAMERA <- 
  function(xset, ...) {
    
    ids <- xset@peaks[,"sample", drop=FALSE]
    ids <- cbind(1:length(ids), ids)
    
    iso_mat <- matrix(0, nrow=0, ncol=2)
    
    xsets <- split(xset, unique(xset@peaks[,"sample"]))
    samples <- unique(xset@peaks[,"sample"])
    for(sample in samples) {
      an <- xsAnnotate(xset, sample=sample)
      isos <- findIsotopes(an, ...)@isoID[,c("mpeak", "isopeak"), drop=FALSE]
      #start_id <- ids[ids[,2]==sample,,drop=FALSE][1,1] - 1
      iso_mat <- rbind(iso_mat, matrix(ids[ids[,2]==sample,1][isos], ncol=2))
    }
    
    iso_mat
  }


calcMaximumCarbon <- 
  function(masses) {  
    
    carbon = 12.0
    hydrogen  = 1.0078250170
    CH3 = carbon + 3 * hydrogen
    CH2 = carbon + 2 * hydrogen  
    
    maximum_carbon <- floor((masses-2*CH3)/CH2) + 2
    
  }    

calculateXcmsSet <- function(example_sample, xcmsSetParameters, scanrange, task=1, nSlaves=1) {
  xset <- NULL  
  
  if(is.null(xcmsSetParameters$step)) {     #centWave    
    xset <- xcmsSet(files=example_sample, method="centWave", 
                    peakwidth=c(xcmsSetParameters$min_peakwidth[task], xcmsSetParameters$max_peakwidth[task]),
                    ppm=xcmsSetParameters$ppm[task], noise=xcmsSetParameters$noise[task], 
                    snthresh=xcmsSetParameters$snthresh[task], mzdiff=xcmsSetParameters$mzdiff[task],
                    prefilter=c(xcmsSetParameters$prefilter[task], xcmsSetParameters$value_of_prefilter[task]),
                    mzCenterFun=xcmsSetParameters$mzCenterFun[task], integrate=xcmsSetParameters$integrate[task],
                    fitgauss=xcmsSetParameters$fitgauss[task], verbose.columns=xcmsSetParameters$verbose.columns[task],
                    scanrange=scanrange, nSlaves=nSlaves)
    
  } else {     #matchedFilter  
    try(xset <- xcmsSet(files=example_sample, method="matchedFilter", 
                        fwhm=xcmsSetParameters$fwhm[task], snthresh=xcmsSetParameters$snthresh[task],
                        step=xcmsSetParameters$step[task], steps=xcmsSetParameters$steps[task],
                        sigma=xcmsSetParameters$sigma[task], max=xcmsSetParameters$max[task], 
                        mzdiff=xcmsSetParameters$mzdiff[task], index=xcmsSetParameters$index[task],
                        scanrange=scanrange, nSlaves=nSlaves))   
  }
  
  return(xset)
  
}


checkXcmsSetParams <-
  function(params) {
    if(is.null(params$step)) {     #centWave   
      quantitative_parameters <- c("ppm", "min_peakwidth", "max_peakwidth", "snthresh", "mzdiff", "noise", "prefilter", "value_of_prefilter")
      qualitative_parameters <- c("integrate", "fitgauss", "verbose.columns", "mzCenterFun")
      unsupported_parameters <- c("sleep", "ROI.list") #"scanrange" can only be set, but not optimized
    } else {    
      quantitative_parameters <- c("fwhm", "sigma", "max", "snthresh", "step", "steps", "mzdiff")
      qualitative_parameters <- c("index")
      unsupported_parameters <- c("sleep")  
    } 
    checkParams(params, quantitative_parameters, qualitative_parameters, unsupported_parameters)
  }


getDefaultXcmsSetStartingParams <-
  function(method=c("centWave", "matchedFilter")) {
    
    method <- match.arg(method)
    
    if(method=="centWave")
      return(list(min_peakwidth=c(12,28), max_peakwidth=c(35,65), ppm=c(17,32),
                  mzdiff=c(-0.001, 0.01), snthresh=10, noise=0, prefilter=3, 
                  value_of_prefilter=100,  mzCenterFun="wMean", integrate=1, 
                  fitgauss=FALSE, verbose.columns=FALSE))
    
    if(method=="matchedFilter")
      return(list(fwhm=c(25,35), snthresh=c(3,17), step=c(0.05, 0.15), steps=c(1,3), 
                  sigma=0, max=5, mzdiff=0, index=FALSE)) 
    
  }


optimizeSlaveCluster <-
  function(task, xcmsSet_parameters, files, scanrange, isotopeIdentification, ...) {
    
    print(task)  
    
    #library(xcms)
    xset <- NULL
    print(sapply(xcmsSet_parameters, "[[", task))
    
    
    xset <- calculateXcmsSet(files, xcmsSet_parameters, scanrange, task)
    
#     if(is.null(xcmsSet_parameters$step)) {     #centWave  	
#       xset <- xcmsSet(files=example_sample, method="centWave", 
#                       peakwidth=c(xcmsSet_parameters$min_peakwidth[task], xcmsSet_parameters$max_peakwidth[task]),
#                       ppm=xcmsSet_parameters$ppm[task], noise=xcmsSet_parameters$noise[task], 
#                       snthresh=xcmsSet_parameters$snthresh[task], mzdiff=xcmsSet_parameters$mzdiff[task],
#                       prefilter=c(xcmsSet_parameters$prefilter[task], xcmsSet_parameters$value_of_prefilter[task]),
#                       mzCenterFun=xcmsSet_parameters$mzCenterFun[task], integrate=xcmsSet_parameters$integrate[task],
#                       fitgauss=xcmsSet_parameters$fitgauss[task], verbose.columns=xcmsSet_parameters$verbose.columns[task])
#       
#     } else {     #matchedFilter  
#       try(xset <- xcmsSet(files=example_sample, method="matchedFilter", 
#                           fwhm=xcmsSet_parameters$fwhm[task], snthresh=xcmsSet_parameters$snthresh[task],
#                           step=xcmsSet_parameters$step[task], steps=xcmsSet_parameters$steps[task],
#                           sigma=xcmsSet_parameters$sigma[task], max=xcmsSet_parameters$max[task], 
#                           mzdiff=xcmsSet_parameters$mzdiff[task], index=xcmsSet_parameters$index[task]))   
#     }
    result <- calcPPS(xset, isotopeIdentification, ...)
    result[1] <- task   
    rm(xset)
    
    result
    
  }


optimizeXcmsSet <-
  function(files=NULL, params=getDefaultXcmsSetStartingParams(), isotopeIdentification=
             c("IPO", "CAMERA"), nSlaves=4, subdir="IPO", ...) { #ppm=5, rt_diff=0.02, nSlaves=4, subdir="IPO") {
    
    scanrange <- params$scanrange
    params$scanrange <- NULL    
    checkXcmsSetParams(params)

    
    if(is.null(files)) {
      files <- getwd()
    }
    
    isotopeIdentification <- match.arg(isotopeIdentification)
    
    centWave <- is.null(params$fwhm)  
    
    history <- list()
    iterator = 1 
    best_range <- 0.25
    
    if(!is.null(subdir))
      if(!file.exists(subdir))
        dir.create(subdir)
    
    
    while(iterator < 50) { #dummy stop :-)
      cat("\n")
      cat("\n")
      cat("\n")
      cat("starting new DoE with:\n")
      print(params)
      
      #    xcms_result <- xcmsSetExperiments(files, params, nSlaves) 
      #                       ppm, rt_diff, nSlaves)   
      xcms_result <- xcmsSetExperimentsCluster(files, params, scanrange, isotopeIdentification, nSlaves, ...) 
      
      
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
        
#         if(centWave) {		
#           xset <- xcmsSet(files=files, method="centWave", 
#                           peakwidth=c(xcms_parameters$min_peakwidth, xcms_parameters$max_peakwidth),
#                           ppm=xcms_parameters$ppm, noise=xcms_parameters$noise, 
#                           snthresh=xcms_parameters$snthresh, mzdiff=xcms_parameters$mzdiff,
#                           prefilter=c(xcms_parameters$prefilter, xcms_parameters$value_of_prefilter),
#                           mzCenterFun=xcms_parameters$mzCenterFun, integrate=xcms_parameters$integrate,
#                           fitgauss=xcms_parameters$fitgauss, verbose.columns=xcms_parameters$verbose.columns, 
#                           scanrange=xcms_parameters$scanrange, nSlaves=nSlaves)
#         } else {
#           xset <- xcmsSet(files=files, method="matchedFilter", 
#                           fwhm=xcms_parameters$fwhm, snthresh=xcms_parameters$snthresh,
#                           step=xcms_parameters$step, steps=xcms_parameters$steps,
#                           sigma=xcms_parameters$sigma, max=xcms_parameters$max, 
#                           mzdiff=xcms_parameters$mzdiff, index=xcms_parameters$index,
#                           scanrange=xcms_parameters$scanrange, nSlaves=nSlaves) 
#           
#           
#         }
        xset <- calculateXcmsSet(files, xcms_parameters, scanrange, 1, nSlaves)
        
        best_settings$xset <- xset
        target_value <- calcPPS(xset, isotopeIdentification, ...)#, ppm, rt_diff)
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
        
        if(names(params$to_optimize)[i] == "steps" | names(params$to_optimize)[i] == "prefilter") {
          params$to_optimize[[i]] <- round(new_bounds, 0)
        } else { 
          params$to_optimize[[i]] <- new_bounds
        }
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
      
      params <- attachList(params$to_optimize, params$no_optimization)	    
      iterator <- iterator + 1
      
    }
    params <- attachList(params$to_optimize, params$no_optimization)	    
    return(history)
    
  }


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



xcmsSetExperimentsCluster <-
function(example_sample, params, scanrange, isotopeIdentification, nSlaves=4, ...) { 

  typ_params <- typeCastParams(params) 

  if(length(typ_params$to_optimize)>1) {
    design <- getCcdParameter(typ_params$to_optimize)  	
    xcms_design <- decode.data(design) 
  } else {
    design <- data.frame(run.order=1:9, a=seq(-1,1,0.25))
      colnames(design)[2] <- names(typ_params$to_optimize)
      xcms_design <- design
      xcms_design[,2] <- seq(min(typ_params$to_optimize[[1]]), max(typ_params$to_optimize[[1]]), diff(typ_params$to_optimize[[1]])/8)
  }

  xcms_design <- combineParams(xcms_design, typ_params$no_optimization)   
  tasks <- 1:nrow(design)  
  
  if(nSlaves > 1) {
    cl <- makeCluster(nSlaves, type = "PSOCK")
    response <- matrix(0, nrow=length(design[[1]]), ncol=5)
 
    #exporting all functions to cluster but only calcPPS and toMatrix are needed
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    clusterExport(cl, ex)
    response <- parSapply(cl, tasks, optimizeSlaveCluster, xcms_design, example_sample, 
                          scanrange, isotopeIdentification, ..., USE.NAMES=FALSE)
    stopCluster(cl)
  } else {
    response <- sapply(tasks, optimizeSlaveCluster, xcms_design, example_sample, 
                       scanrange, isotopeIdentification, ...)
  }
  
  response <- t(response)
  colnames(response) <- c("exp", "num_peaks", "notLLOQP", "num_C13", "PPS")
  response <- response[order(response[,1]),]

  ret <- list()
  ret$params <- typ_params
  ret$design <- design
  ret$response <- response

  return(ret)

}


xcmsSetStatistic <-
function(xcms_result, subdir, iterator) {

  params <- xcms_result$params
  resp <- xcms_result$response[, "PPS"]
  
 
    model <- createModel(xcms_result$design, params$to_optimize, resp)
    xcms_result$model <- model                  
     
    max_settings <- getMaximumLevels(xcms_result$model)
    tmp <- max_settings[1,-1]
    tmp[is.na(tmp)] <- 1
    if(!is.null(subdir) & length(tmp) > 1)
      plotContours(xcms_result$model, tmp, paste(subdir,"/rsm_", iterator, sep=""))
	
  xcms_result$max_settings <- max_settings

  return(xcms_result)
}

