startSlaves <-
function(n_slaves) {
  mpi.spawn.Rslaves(nslaves=n_slaves)
                                                                                
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
