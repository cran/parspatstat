# Creation of parppp class
#  nrow - number of rows in full grid
#  ncol - number of cols in full grid
#  ranks - vector of ranks corresponding to data
#  topology - matrix ilustrating how data should be combined
#  ncpu - number of CPUs (redundant?)
#  data - list of SpatialPixelsDataFrame or SpatialGridDataFrame in rank order
setClass("parppp",
         representation(
           localname="character", 
           nslaves="numeric",
           comm="numeric"),
         prototype = list(
           localname="X.ppp",
           nslaves=numeric(),
           comm=numeric())
)


print.parppp <- function(x, ...) {
  print(paste("Parallel point pattern object named",x@localname,"in",x@nslaves,"slaves on comm",x@comm))
  print("Point pattern breakdown:")
  ppps <- mpi.remote.exec("get",x@localname,envir=.BaseNamespaceEnv,comm=x@comm)
  subwindows <- mpi.remote.exec("get",paste(x@localname,".subwin",sep=""),envir=.BaseNamespaceEnv,comm=x@comm)
  for (i in 1:x@nslaves) {
    print(paste("Slave",i,"-",ppps[[i]]$n,"points"))
    print(subwindows[[i]])
  }
  print("Full window")
  print(ppps[[1]]$window)
}


plot.parppp <- function(x, ...) {
  ppps <- mpi.remote.exec("get",x@localname,envir=.BaseNamespaceEnv,comm=x@comm)

  for (i in 1:x@nslaves) {
    if (i==1) add=FALSE
    else add=TRUE
    plot(ppps[[i]],add=add,...)
  }
  
}