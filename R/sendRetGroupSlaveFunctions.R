sendRetGroupSlaveFunctions <-
function(parameters, xset) {
  mpi.bcast.Robj2slave(getRTGVValues)
  mpi.bcast.cmd(slave <- mpi.comm.rank())
  mpi.bcast.Robj2slave(parameters)
  mpi.bcast.Robj2slave(xset)  
  mpi.bcast.cmd("library(xcms)")
  mpi.bcast.Robj2slave(optimizeRetGroupSlave)  
  mpi.bcast.cmd(optimizeRetGroupSlave())
}
