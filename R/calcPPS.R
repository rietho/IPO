calcPPS <-
function(xset) {
  
  peak_source <- xset@peaks[,c("mz", "rt", "sample", "into", "mzmin", "mzmax", "rtmin", "rtmax")]
   
  iso_list <- list()
  ret <- array(0, dim=c(1,5))  
  if(nrow(peak_source) == 0) {
	  return(ret)
  }  
  ret[2] <- nrow(xset@peaks)

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
	  speaks <- peak_source[peak_source[,"sample"]==sample,]	
	  found_isotope <- FALSE
    split <- 250
    	
	  if(!is.null(ncol(speaks)) & length(speaks) > 3) {  		      
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
	        part_peaks <- speaks[1:end_point,]
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

      sample_isos_peaks <- xset@peaks
      sample_non_isos_peaks <- xset@peaks
      
      if(found_isotope) {
        sample_isos_peaks <- xset@peaks[unique(unlist(iso_list)),]
        sample_non_isos_peaks <- xset@peaks[-unique(unlist(iso_list)),]
      } 

      speaks <- sample_non_isos_peaks[sample_non_isos_peaks[,"sample"]==sample,]
      sample_isos_peaks <- sample_isos_peaks[sample_isos_peaks[,"sample"]==sample,]
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
      ret[5] <- ret[4]^1.5/ret[3]  
	  }    
  }#end_for_sample    
  
  return(ret)

}
