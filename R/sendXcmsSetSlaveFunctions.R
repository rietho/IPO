sendXcmsSetSlaveFunctions <-
function(example_sample, xcmsSet_parameters) {#, ppm, rt_diff) {
  mpi.bcast.Robj2slave(toMatrix) 
  mpi.bcast.Robj2slave(calcPPS) 
  mpi.bcast.cmd(slave <- mpi.comm.rank())
  mpi.bcast.Robj2slave(example_sample)
  mpi.bcast.Robj2slave(xcmsSet_parameters)
  mpi.bcast.Robj2slave(optimizeSlave)
  
  mpi.bcast.cmd("library(xcms)")
  mpi.bcast.cmd(optimizeSlave())
}
