plotContours <-
function(model, maximum_slice, plot_name) {

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
