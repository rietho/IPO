calcPPS <- function(xset, isotopeIdentification=c("IPO", "CAMERA"), ...) {
    
    isotopeIdentification <- match.arg(isotopeIdentification)
    
    ret <- vector(mode="numeric", 5) #array(0, dim=c(1,5)) 
    names(ret) <- c("ExpId", "#peaks", "#NonRP", "#RP", "PPS")
    if(is.null(xset)) {
      return(ret)
    } 
    
    if(nrow(peaks(xset)) == 0) {
      return(ret)
    }
    
    peak_source <- peaks(xset)[,c("mz", "rt", "sample", "into", "mzmin", 
                                 "mzmax", "rtmin", "rtmax"),drop=FALSE]
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
        non_isos_peaks <- peak_source[-unique(c(iso_mat)),,drop=FALSE] 
      } 
	  
      speaks <- non_isos_peaks[non_isos_peaks[,"sample"]==sample,,drop=FALSE]
      intensities <- speaks[,"into"]
	    na_int <- is.na(intensities)
	    intensities <- intensities[!na_int]
      
	    if(length(intensities)>0) {
        tmp <- intensities[order(intensities)]
        int_cutoff <- mean(tmp[1:max(round((length(tmp)/33),0),1)])
      
        masses <- speaks[!na_int, "mz"]
        #floor((masses-2*CH3)/CH2) + 2
        maximum_carbon <- calcMaximumCarbon(masses)
        carbon_probabilty <- maximum_carbon*isotope_abundance
      
        iso_int <- intensities * carbon_probabilty
      
        not_loq_peaks <- sum(iso_int>int_cutoff)
        ret[3] <- ret[3] + not_loq_peaks
      }
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
  function(xset, checkPeakShape=c("none", "borderIntensity", "sinusCurve", 
                                  "normalDistr")) {
    
    checkPeakShape <- match.arg(checkPeakShape)
    
    iso_mat <- matrix(0, nrow=0, ncol=2)
    if(is.null(xset)) {
      return(iso_mat)
    }
    
    colnames(iso_mat) <- c("12C", "13C")
    peak_source <- peaks(xset)[,c("mz", "rt", "sample", "into", "maxo", "mzmin",
                                 "mzmax", "rtmin", "rtmax"), drop=FALSE]
    
    for(i in 1:ncol(peak_source)) {
      peak_source <- peak_source[!is.na(peak_source[,i]),,drop=FALSE]
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
      speaks <- peak_source[peak_source[,"sample"]==sample,,drop=FALSE]
      split <- 250
	  if(!(checkPeakShape=="none"))
        rawdata <- loadRaw(xcmsSource(filepaths(xset)[sample]))
      
      if(nrow(speaks)>1) {  		      
        #speaks <- speaks[,-c("sample")]
        speaks <- speaks[order(speaks[,"mz"]),]
        
        while(!is.null(nrow(speaks)) & length(speaks) > 3) {
          part_peaks <- NULL
          #splitting the data into smaller pieces to improve speed    
          if(nrow(speaks) < split) {
            part_peaks <- speaks
          } else {          
            upper_bound <- speaks[split,"mzmax"] + isotope_mass          
            end_point <- sum(speaks[,"mz"] < upper_bound)
            part_peaks <- speaks[1:end_point,,drop=FALSE]
          }		
          
          rt <- part_peaks[,"rt"]
          rt_window <- rt * 0.005
          rt_lower <- part_peaks[,"rt"] - rt_window
          rt_upper <- part_peaks[,"rt"] + rt_window
          rt_matrix <-  
            t(matrix(rep(rt, nrow(part_peaks)), ncol=nrow(part_peaks)))
          rt_matrix_bool <- rt_matrix >= rt_lower & rt_matrix <= rt_upper
          
          mz <- part_peaks[,"mz"]
          #isotope_masses - mz_window
          mz_lower <- part_peaks[,"mzmin"] + isotope_mass
          #isotope_masses + mz_window
          mz_upper <- part_peaks[,"mzmax"] + isotope_mass
          mz_matrix <-  
            t(matrix(rep(mz, nrow(part_peaks)), ncol=nrow(part_peaks)))
          mz_matrix_bool <- mz_matrix >= mz_lower & mz_matrix <= mz_upper
          
          rt_mz_matrix_bool <- rt_matrix_bool & mz_matrix_bool
          
          rt_mz_peak_ids <- which(rowSums(rt_mz_matrix_bool)>0)
          calculations <- min(split, nrow(speaks))
          rt_mz_peak_ids <- rt_mz_peak_ids[rt_mz_peak_ids < calculations]
          
          for(i in rt_mz_peak_ids) {
            current <- part_peaks[i, ,drop=FALSE]
            rt_mz_peaks <- part_peaks[rt_mz_matrix_bool[i,],,drop=FALSE]
            rt_difference <- 
              abs(current[,"rt"] - rt_mz_peaks[, "rt"]) / current[,"rt"]
            rt_mz_peaks <- cbind(rt_mz_peaks, rt_difference)
            #test intensity_window
            #floor((current["mz"]-2*CH3)/CH2) + 2
            maximum_carbon <- calcMaximumCarbon(current[,"mz"]) 
            carbon_probabilty <- c(1,maximum_carbon)*isotope_abundance
            iso_intensity <- current[,"into"] * carbon_probabilty
            
            int_bools <- 
              rt_mz_peaks[,"into"] >= iso_intensity[1] & 
              rt_mz_peaks[,"into"] <= iso_intensity[2]
            
            if(sum(int_bools) > 0) {
              int_peaks <- rt_mz_peaks[int_bools,,drop=FALSE]
              boundary_bool <- rep(TRUE, (nrow(int_peaks)+1))
              if(!(checkPeakShape=="none")) {
                if(checkPeakShape=="borderIntensity") {
                  boundary_bool <- checkIntensitiesAtRtBoundaries(
                    rawdata, 
                    rbind(current,int_peaks[,-ncol(int_peaks), drop=FALSE]))
                } else {
                  if(checkPeakShape=="sinusCurve") {                
                    boundary_bool <- checkSinusDistribution(
                      rawdata, 
                      rbind(current,int_peaks[,-ncol(int_peaks),drop=FALSE]))
                  } else {                  
                    boundary_bool <- checkNormalDistribution(
                      rawdata, 
                      rbind(current,int_peaks[,-ncol(int_peaks),drop=FALSE]))
                  }
                }
              } #else {
                #boundary_bool <- rep(TRUE, (nrow(int_peaks)+1))
              #}              
              if(boundary_bool[1] & sum(boundary_bool[-1])>0) {                 
                iso_peaks <- int_peaks[boundary_bool[-1],,drop=FALSE]
                iso_id <- 
                  iso_peaks[which.min(iso_peaks[,"rt_difference"]), "id"]
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

# checking intensities at rtmin and rtmax. peaks[,"maxo"] must be at least 
# double as high does not work for retention time corrected data 
checkIntensitiesAtRtBoundaries <-
  function(rawdata, 
           peaks, 
           minBoundaryToMaxo=1/3, 
           ppmstep=15) {
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
        mzIntensities <- 
          c(0, intensities[mz>=peak["mzmin"]-ppm[1] & mz<=peak["mzmax"]+ppm[2]])
        maxBoundaryIntensity <- max(mzIntensities)
        ret[i] <- ret[i] & maxBoundaryIntensity<peak["maxo"]*minBoundaryToMaxo
      }
    }
  }
  
  return(ret)
  
}

checkSinusDistribution <- function(rawdata, peaks) {
 ret <- rep(TRUE, nrow(peaks))
 for(i in 1:nrow(peaks)) {
   ret[i] <- testSinusDistribution(rawdata, peaks[i,,drop=FALSE])
 }
 
 return(ret)
}

checkNormalDistribution <- function(rawdata, peaks) {
 ret <- rep(TRUE, nrow(peaks))
 for(i in 1:nrow(peaks)) {
   ret[i] <- testNormalDistribution(rawdata, peaks[i,,drop=FALSE])
 }
 
 return(ret)
}

getIntensitiesFromRawdata <- function(rawdata, peak) {
  rt <- rawdata$rt >= peak[,"rtmin"] & rawdata$rt <= peak[,"rtmax"]

  rtRange <- c(min(which(rt)), max(which(rt))+1)  
  scanIndices <- 
    rawdata$scanindex[rtRange[1]:min(rtRange[2], length(rawdata$scanindex))]
  #  scanIndices <- scanIndices[!is.na(scanIndices)]
  if(rtRange[2]>length(rawdata$scanindex)) {
    scanIndices <- c(scanIndices, length(rawdata$intensity))
  }
  
  if(length(scanIndices) < 3)
    return(FALSE)  
  
  y <- c()
  for(i in 1:(length(scanIndices)-1)) {
    scanRange <- c(scanIndices[i]+1, scanIndices[i+1])
    mz <- rawdata$mz[scanRange[1]:scanRange[2]]
    y <- 
      c(y, 
        max(0, (rawdata$intensity[scanRange[1]:scanRange[2]][
          mz >= peak[,"mzmin"] & mz <= peak[,"mzmax"]])
            )
        )
  }
  
  y
}

testNormalDistribution <- function(rawdata, peak) {

  y <- getIntensitiesFromRawdata(rawdata, peak)
  if(length(y) < 3) {
    return(FALSE)
  }
  
  if(max(y)==0) {
    return(FALSE)
  }
  
  normY <- (y-min(y))/(max(y)-min(y))
  
  mean=10; 
  sd=3;

  seqModel <- seq(-4,4,length=length(normY))*sd + mean
  yModel <- dnorm(seqModel,mean,sd)
  yModel = yModel* (1/max(yModel))
  correlation <- cor(yModel, normY)

  correlation > 0.7

  
}

testSinusDistribution <- function(rawdata, peak) {

  y <- getIntensitiesFromRawdata(rawdata, peak)
  if(length(y) < 3) {
    return(FALSE)
  }
  if(max(y)==0) {
    return(FALSE)
  }

  normY <- (y-min(y))/(max(y)-min(y))
  sinCurve <- (sin(seq(-pi/2,pi+1.5,length=length(normY))) + 1) / 2
  correlation <- cor(sinCurve, normY)

  correlation > 0.7
  
}

findIsotopes.CAMERA <- 
  function(xset, ...) {
    
    iso_mat <- matrix(0, nrow=0, ncol=2)
    if(is.null(xset)) {
      return(iso_mat)
    }
    
    ids <- peaks(xset)[,"sample", drop=FALSE]
    ids <- cbind(1:length(ids), ids)
    
    xsets <- split(xset, unique(peaks(xset)[,"sample"]))
    samples <- unique(peaks(xset)[,"sample"])
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

calculateXcmsSet <- 
  function(files, xcmsSetParameters, scanrange=NULL, task=1, nSlaves=1) {
  xset <- NULL  
  
  if(is.null(xcmsSetParameters$step)) { # centWave    
    xset <- 
      xcmsSet(files = files, 
              method = "centWave", 
              peakwidth = 
                c(xcmsSetParameters$min_peakwidth[task],
                  xcmsSetParameters$max_peakwidth[task]),
              ppm         = xcmsSetParameters$ppm[task], 
              noise       = xcmsSetParameters$noise[task], 
              snthresh    = xcmsSetParameters$snthresh[task], 
              mzdiff      = xcmsSetParameters$mzdiff[task],
              prefilter = 
                c(xcmsSetParameters$prefilter[task],
                  xcmsSetParameters$value_of_prefilter[task]),
              mzCenterFun = xcmsSetParameters$mzCenterFun[task], 
              integrate   = xcmsSetParameters$integrate[task],
              fitgauss    = xcmsSetParameters$fitgauss[task], 
              verbose.columns = 
                xcmsSetParameters$verbose.columns[task],
              scanrange   = scanrange, 
              nSlaves     = nSlaves*xcmsSetParameters$nSlaves[task])
    
  } else {     #matchedFilter  
    try(xset <- 
          xcmsSet(files     = files, 
                  method    = "matchedFilter", 
                  fwhm      = xcmsSetParameters$fwhm[task], 
                  snthresh  = xcmsSetParameters$snthresh[task],
                  step      = xcmsSetParameters$step[task], 
                  steps     = xcmsSetParameters$steps[task],
                  sigma     = xcmsSetParameters$sigma[task], 
                  max       = xcmsSetParameters$max[task], 
                  mzdiff    = xcmsSetParameters$mzdiff[task], 
                  index     = xcmsSetParameters$index[task],
                  scanrange = scanrange, 
                  nSlaves   = nSlaves*xcmsSetParameters$nSlaves[task]))  
  }
  return(xset)
}


checkXcmsSetParams <-
  function(params) {
    if(is.null(params$step)) { # centWave   
      quantitative_parameters <- 
        c("ppm", "min_peakwidth", "max_peakwidth", "snthresh", 
          "mzdiff", "noise", "prefilter", "value_of_prefilter")
      qualitative_parameters <- 
        c("integrate", "fitgauss", "verbose.columns", "mzCenterFun",
          "nSlaves")
      unsupported_parameters <- 
        c("sleep", "ROI.list") # "scanrange" can only be set,
                               # but not optimized
    } else {    
      quantitative_parameters <- 
        c("fwhm", "sigma", "max", "snthresh", "step", "steps", "mzdiff")
      qualitative_parameters <- 
        c("index", "nSlaves")
      unsupported_parameters <- 
        c("sleep")  
    } 
    checkParams(params, quantitative_parameters, qualitative_parameters,
                unsupported_parameters)
  }


getDefaultXcmsSetStartingParams <-
  function(method=c("centWave", "matchedFilter")) {
    
    method <- match.arg(method)
    
    if(method=="centWave")
      return(list(min_peakwidth      = c(12,28), 
                  max_peakwidth      = c(35,65), 
                  ppm                = c(17,32),
                  mzdiff             = c(-0.001, 0.01), 
                  snthresh           = 10, 
                  noise              = 0, 
                  prefilter          = 3, 
                  value_of_prefilter = 100,  
                  mzCenterFun        = "wMean", 
                  integrate          = 1, 
                  fitgauss           = FALSE, 
                  verbose.columns    = FALSE, 
                  nSlaves            = 1))
    
    if(method=="matchedFilter")
      return(list(fwhm     = c(25,35), 
                  snthresh = c(3,17), 
                  step     = c(0.05, 0.15), 
                  steps    = c(1,3), 
                  sigma    = 0,
                  max      = 5, 
                  mzdiff   = 0, 
                  index    = FALSE, 
                  nSlaves  = 1)) 
    
  }


optimizeSlaveCluster <-
  function(task, 
           xcmsSet_parameters,
           files, 
           scanrange, 
           isotopeIdentification, 
           ...) {
    
    message(task)  
    
    #library(xcms)
    xset <- NULL
    #message(sapply(xcmsSet_parameters, "[[", task))
    
    
    xset <- calculateXcmsSet(files, xcmsSet_parameters, scanrange, task)
    
    result <- calcPPS(xset, isotopeIdentification, ...)
    result[1] <- task   
    rm(xset)
    
    result
    
  }


optimizeXcmsSet <-
  function(files=NULL, params=getDefaultXcmsSetStartingParams(), 
           isotopeIdentification=c("IPO", "CAMERA"), nSlaves=4, 
           subdir="IPO", ...) { 
    
    scanrange <- params$scanrange
    params$scanrange <- NULL    
    checkXcmsSetParams(params)

    
    if(is.null(files)) {
      files <- getwd()
    }
    
    isotopeIdentification <- match.arg(isotopeIdentification)
    
    # only matchedFilter-method has parameter fwhm
    centWave <- is.null(params$fwhm)  
    
    history <- list()
    iterator = 1 
    best_range <- 0.25
    
    if(!is.null(subdir))
      if(!file.exists(subdir))
        dir.create(subdir)
    
    
    while(iterator < 50) { #dummy stop :-)
      message("\n\n")
      message("starting new DoE with:")
      message(paste(rbind(paste(names(params), sep="", ": "), 
                          paste(params, sep="", "\n")),sep=""))
      
      xcms_result <- 
        xcmsSetExperimentsCluster(files, params, scanrange, 
                                  isotopeIdentification, nSlaves, ...) 
	  
      xcms_result <- xcmsSetStatistic(files, scanrange, 
                                      isotopeIdentification, xcms_result, 
                                      subdir, iterator, nSlaves, ...)
      history[[iterator]] <- xcms_result     
      params <- xcms_result$params 
      
      if(!resultIncreased(history)) {
        message("no increase, stopping")
        maxima <- 0
        max_index <- 1
        for(i in 1:length(history)) {
          if(history[[i]]$max_settings[1] > maxima) {
            maxima <- history[[i]]$max_settings[1]
            max_index <- i
          }
        }
        
        xcms_parameters <- 
          as.list(decodeAll(history[[max_index]]$max_settings[-1],
                            history[[max_index]]$params$to_optimize))      
        xcms_parameters <- combineParams(xcms_parameters, 
                                         params$no_optimization)
        
        if(!is.list(xcms_parameters))
          xcms_parameters <- as.list(xcms_parameters)
        
        best_settings <- list()
        best_settings$parameters <- xcms_parameters
        
        best_settings$xset <- history[[max_index]]$xset
        #calcPPS(xset, isotopeIdentification, ...)#, ppm, rt_diff)
        target_value <- history[[max_index]]$PPS 
        best_settings$result <- target_value
        history$best_settings <- best_settings
        
        message("best parameter settings:")
        message(paste(rbind(paste(names(xcms_parameters), 
                                  sep="", ": "), 
                            paste(xcms_parameters, sep="", "\n")), sep=""))
        
        return(history)   
      }
      
      
      for(i in 1:length(params$to_optimize)) {
        parameter_setting <- xcms_result$max_settings[i+1]
        bounds <- params$to_optimize[[i]] 
        fact <- names(params$to_optimize)[i]
        min_factor <- 
          ifelse(fact=="min_peakwidth", 3, 
                 ifelse(fact=="mzdiff", 
                        ifelse(centWave,-100000000, 0.001), 
                             ifelse(fact=="step",0.0005,1)))
        
        # if the parameter is NA, we increase the range by 20%, 
        # if it was within the inner 
        # 25% of the previous range or at the minimum value we decrease 
        # the range by 20%
        step_factor <- 
          ifelse(is.na(parameter_setting), 1.2, 
                 ifelse((abs(parameter_setting) < best_range), 
                              0.8, 
                        ifelse(parameter_setting==-1 & 
                                 decode(-1, params$to_optimize[[i]]) ==
                              min_factor, 0.8, 1)))
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
        
        if(names(params$to_optimize)[i] == "steps" | 
           names(params$to_optimize)[i] == "prefilter") {
          params$to_optimize[[i]] <- round(new_bounds, 0)
        } else { 
          params$to_optimize[[i]] <- new_bounds
        }
      } 
      
      if(centWave) {
        #checking peakwidths plausiability
        if(!is.null(params$to_optimize$min_peakwidth) | 
           !is.null(params$to_optimize$max_peakwidth)) {
          pw_min <- 
            ifelse(is.null(params$to_optimize$min_peakwidth), 
                   params$no_optimization$min_peakwidth, 
                   max(params$to_optimize$min_peakwidth))
          pw_max <- 
            ifelse(is.null(params$to_optimize$max_peakwidth), 
                   params$no_optimization$max_peakwidth, 
                   min(params$to_optimize$max_peakwidth))
          if(pw_min >= pw_max) {
            additional <- abs(pw_min-pw_max) + 1
            if(!is.null(params$to_optimize$max_peakwidth)) {		  
              params$to_optimize$max_peakwidth <- 
                params$to_optimize$max_peakwidth + additional
            } else {
              params$no_optimization$max_peakwidth <- 
                params$no_optimization$max_peakwidth + additional
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
    if(history[[index]]$PPS["PPS"] == 0 & index == 1)
      stop(paste("No isotopes have been detected,",
                 "peak picking not optimizable by IPO!"))
    
    if(index < 2)
      return(TRUE)
    
    if(history[[index-1]]$PPS["PPS"] >= history[[index]]$PPS["PPS"])
      return(FALSE)
    
    return(TRUE)
    
  }



xcmsSetExperimentsCluster <-
  function(example_sample, 
           params, 
           scanrange, 
           isotopeIdentification, 
           nSlaves=4, 
           ...) { 

  typ_params <- typeCastParams(params) 

  if(length(typ_params$to_optimize)>1) {
    design <- getCcdParameter(typ_params$to_optimize)  	
    xcms_design <- decode.data(design) 
  } else {
    design <- data.frame(run.order=1:9, a=seq(-1,1,0.25))
      colnames(design)[2] <- names(typ_params$to_optimize)
      xcms_design <- design
      xcms_design[,2] <- 
        seq(min(typ_params$to_optimize[[1]]), 
            max(typ_params$to_optimize[[1]]), 
            diff(typ_params$to_optimize[[1]])/8)
  }

  xcms_design <- combineParams(xcms_design, typ_params$no_optimization)   
  tasks <- 1:nrow(design)  
  
  if(nSlaves > 1) {
    # unload snow (if loaded) to prevent conflicts with usage of
    # package 'parallel'
    if('snow' %in% rownames(installed.packages()))
      unloadNamespace("snow")
    
    cl_type<-getClusterType()
    cl <- parallel::makeCluster(nSlaves, type = cl_type)
    response <- matrix(0, nrow=length(design[[1]]), ncol=5)
 
    #exporting all functions to cluster but only calcPPS and toMatrix are needed
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    if(identical(cl_type,"PSOCK")) {
      message("Using PSOCK type cluster, this increases memory requirements.")
      message("Reduce number of slaves if you have out of memory errors.")
      message("Exporting variables to cluster...")
      parallel::clusterExport(cl, ex)
    }
    response <- parallel::parSapply(cl, tasks, optimizeSlaveCluster,
                                    xcms_design, example_sample, scanrange, 
                                    isotopeIdentification, ..., USE.NAMES=FALSE)
    parallel::stopCluster(cl)
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
  function(files, 
           scanrange, 
           isotopeIdentification, 
           xcms_result, 
           subdir, 
           iterator,  
           nSlaves, 
           ...) {

  params <- xcms_result$params
  resp <- xcms_result$response[, "PPS"]
  
 
  model <- createModel(xcms_result$design, params$to_optimize, resp)
  xcms_result$model <- model                  
     
  max_settings <- getMaximumLevels(xcms_result$model)
  
  # plotting rms
  tmp <- max_settings[1,-1] # first row without response
  tmp[is.na(tmp)] <- 1 # if Na (i.e. -1, and 1), show +1
  if(!is.null(subdir) & length(tmp) > 1)
    plotContours(xcms_result$model, tmp, 
                 paste(subdir, "/rsm_", iterator, sep=""))
	
  xcms_result$max_settings <- max_settings
  
  xcms_parameters <- 
    as.list(decodeAll(max_settings[-1], params$to_optimize))      
  xcms_parameters <- 
    combineParams(xcms_parameters, params$no_optimization)
        
  if(!is.list(xcms_parameters))
     xcms_parameters <- as.list(xcms_parameters)
  
  xset <- calculateXcmsSet(files, xcms_parameters, scanrange, 1, nSlaves)
  xcms_result$xset <- xset
  xcms_result$PPS <- calcPPS(xset, isotopeIdentification, ...)  

  return(xcms_result)
}

