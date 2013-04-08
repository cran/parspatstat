parKest <- function(X, ...) {
  if (is.null(X)) {
    stop("Error: No data supplied; data must be of type ppp or parppp")
  } 
  
  #apply the appropriate method
  if (class(X) == "ppp") {
    result <- parKest.ppp(X, ...)
  } else if (class(X)=="parppp") {
    result <- parKest.parppp(X, ...)
  } else {
    stop("X is not an appropriate point process pattern (ppp or parppp) object.")
  }
  
  return(result)
}



parKest.ppp <- function(X, ..., job.num=2*(mpi.comm.size(comm)-1), comm=1, verbose=FALSE, load.sort=TRUE) {
#                     ..., r=NULL, breaks=NULL, 
#                      correction=c("border", "isotropic", "Ripley", "translate"),
#                     nlarge=Inf, domain=NULL, var.approx=FALSE, ratio=FALSE) {
  
  # Error checking
  if (!is.null(X)) {
    if (!verifyclass(X,"ppp")) stop("X is not an appropriate point process pattern (ppp) object.")
  } else {
    stop("Error: No data supplied; data must be of type ppp or parppp.")
  }
  
  nslaves <- mpi.comm.size(comm) - 1
	if(verbose) print(paste(nslaves,"slaves found running on comm",comm))
  if (nslaves < 1) stop("Error: There are no slaves running.")
  if (job.num < 2) stop("Error: job.num needs to be at least 2.")
  if (job.num < nslaves) stop("Error: job.num needs equal or larger than number of running slaves")

	# Send the comm to the slaves
	mpi.bcast.Robj2slave(comm)
	mpi.bcast.Robj2slave(verbose)
	
  if(verbose) print("No errors found, sending slaves into loops") 

	# Have all slaves go into a loop to await commands from master
  mpi.bcast.cmd(cmd=slave.parKest.ppp, X=X, job.num=job.num, ..., verbose=verbose, communicator=comm, comm=comm)

  # Send number of jobs to all slaves (used for tag tracking)
#  mpi.bcast(as.integer(job.num), type=1, comm=comm)
  
  # Send the ppp and arguments to each slave
  #mpi.bcast.Robj(list(X=X, dot.arg=list(...)), rank=0, comm=comm)
 
  if(verbose) print("Objects sent to slaves, awaiting OK status")
		
	# Receive an OK status from slaves before continuting
	for (i in 1:nslaves) {
		tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
		srctag <- mpi.get.sourcetag()

		if(verbose) print(paste("--Slave",srctag[1],"ready to receive jobs"))
	}

  if(verbose) print("Global information sent to slaves, dividing subwindows")
  
  # Create a list of all subwindows
  subwin <- owin.divide(w=X$window, n=job.num)
  if(verbose) print("Subwindows divided, sending first jobs to slaves")
  
  # Find out how many points per job if sorting for load balancing
  if (load.sort) {
    numpointper <- numeric(job.num) # vector to store number of points per subwindow
    for (i in 1:job.num) {
      numpointper[i] <- X[subwin[[i]]]$n
    }
    
    # Sort the jobs by size in each window from highest to lowest
    job.order <- order(numpointper,decreasing=TRUE)
    
    if(verbose) print("Subwindows sorted by decreasing difficulty")
  } else {
    job.order <- 1:job.num #otherwise just default job ordering
  }
 
  # Send first nslaves job to slaves
  for (i in 1:nslaves) {
    mpi.send.Robj(subwin[[job.order[i]]], dest=i, tag=job.order[i], comm=comm)
    if(verbose) print(paste("--Job",job.order[i],"sent to slave", i))
  }
  
  # Loop through all jobs
  out <- as.list(integer(job.num)) #preallocate list to store output
  for (i in 1:job.num) {
    tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
    srctag <- mpi.get.sourcetag()
    
    out[[srctag[2]]] <- tmp # save returned object from slave according to tag number
  		
		if(verbose) print(paste("--Job",srctag[2],"received from slave",srctag[1]))
    
    # Check if all jobs have been sent out
    j <- i + nslaves #the next job
    if (j <= job.num) {
      # if not, send the next job to the slave that just returned
      mpi.send.Robj(obj=subwin[[job.order[j]]], dest=srctag[1], tag=job.order[j], comm=comm)
      if(verbose) print(paste("--Job",job.order[j],"sent to slave", srctag[1]))
    } else {
      # otherwise tell the slave it's finished
      mpi.send.Robj(obj=as.integer(0), dest=srctag[1], tag=0, comm=comm)
			if(verbose) print(paste("--Termination command sent to slave", srctag[1]))
    }
  }

  if(verbose) print("All jobs complete and received, formatting output")
  
  # combine the sub K objects into the whole one (modified from pool.rat)
  Kobj <- pool.K(out)
  
  return(Kobj)
}

slave.parKest.ppp <- function(X, job.num, verbose=FALSE, ..., communicator) {
  #if (!exists("verbose")) verbose=FALSE
  if(verbose) print("Slave initialized, awaiting global information") 

  # get total number of jobs from master
  #job.num <- mpi.bcast(integer(1), type=1, rank=0, comm=1)
  
  # get ppp and other arguments from master
  #tmp <- mpi.bcast.Robj(rank=0, comm=1)
  #X <- tmp$X
  #dotarg <- tmp$dot.arg

	# send ok status to master
	mpi.send.Robj(obj=1, dest=0, tag=1, comm=communicator)
  
  if(verbose) print("Global information received, entering loop to receive jobs")

  #continually get jobs until there are no jobs left (indicated by tag)
  repeat {
    subwin <- mpi.recv.Robj(source=0, tag=mpi.any.tag(), comm=1)
    tag <- mpi.get.sourcetag()[2]
    print(tag)
    if (tag == 0)
      break
    
    #K <- try(do.call("Kest", c(list(X=X, domain=subwin, nlarge=Inf, ratio=TRUE), dotarg)))
    K <- Kest(X=X, domain=subwin, nlarge=Inf, ratio=TRUE, ...)
    attr(attr(K,"fmla"),".Environment") <- emptyenv() #remove unefficiently serialized environment
    
    mpi.send.Robj(obj=K, dest=0, tag=tag, comm=1)
  }
  if(verbose) print("Jobs completed, no more jobs, exiting loop")

}

#parKest for parallel (already divided) point process pattern objects (parppp)
parKest.parppp <- function(X, ..., verbose=FALSE) {
  #                     ..., r=NULL, breaks=NULL, 
  #                      correction=c("border", "isotropic", "Ripley", "translate"),
  #                     nlarge=Inf, domain=NULL, var.approx=FALSE, ratio=FALSE) {
  # Error checking
  if (!is.null(X)) {
    if (!verifyclass(X,"parppp")) stop("X is not an appropriate parallel point process pattern (parppp) object.")
  } else {
    stop("Error: No data supplied; data must be of type parppp.")
  }
  
  nslaves = X@nslaves
  job.num = nslaves
  comm = X@comm
  localname = X@localname
  if(verbose) print(paste(nslaves,"slaves running on comm",comm))
  if (nslaves < 1) stop("Error: There are no slaves running.")
  if (job.num < 2) stop("Error: job.num needs to be at least 2.")
  if (job.num != nslaves) stop("Error: job.num needs to be equal to the number of running slaves")
  
  if(verbose) print("No errors found, sending slaves into function") 
  
  # Have all slaves go into function to perform calcuations from master
  mpi.bcast.cmd(cmd=slave.parKest.parppp, verbose=verbose, localname=localname, ..., communicator=comm, comm=comm)
  
  # Receive an OK status from slaves before continuting
  for (i in 1:nslaves) {
    tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
    srctag <- mpi.get.sourcetag()
    
    if(verbose) print(paste("--Slave",srctag[1],"ready to run"))
  }
  
  job.order <- 1:job.num #just default job ordering since there are only nslaves jobs
  
  # Loop through all jobs
  out <- as.list(integer(job.num)) #preallocate list to store output
  for (i in 1:job.num) {
    #receive K object
    tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
    srctag <- mpi.get.sourcetag()
    
    out[[srctag[2]]] <- tmp # save returned object from slave according to tag number
    
    if(verbose) print(paste("--Job",srctag[2],"received from slave",srctag[1]))
  }
  
  if(verbose) print("All jobs complete and received, formatting output")
  
  # combine the sub K objects into the whole one (modified from pool.rat)
  Kobj <- pool.K(out)
  
  return(Kobj)
}

slave.parKest.parppp <- function(verbose=FALSE, localname, communicator, ...) {
  if(verbose) print("Slave initialized, awaiting global information") 
  
  # send ok status to master
  mpi.send.Robj(obj=1, dest=0, tag=1, comm=communicator)
  
  if(verbose) print("Global information received, entering loop to receive jobs")
  
  #K <- try(do.call("Kest", c(list(X=X, domain=subwin, nlarge=Inf, ratio=TRUE), ...)))
  K <- Kest(X=get(localname,envir=.BaseNamespaceEnv), domain=get(paste(localname,".subwin",sep=""),envir=.BaseNamespaceEnv), nlarge=Inf, ratio=TRUE, ...)
  attr(attr(K,"fmla"),".Environment") <- emptyenv() #remove unefficiently serialized environment
    
  mpi.send.Robj(obj=K, dest=0, tag=mpi.comm.rank(), comm=communicator)
  
  if(verbose) print("Jobs completed, no more jobs, exiting loop")
  
}
