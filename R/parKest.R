parKest <- function(X, ..., job.num=2*(mpi.comm.size(comm)-1), comm=1, verbose=FALSE) {
#                     ..., r=NULL, breaks=NULL, 
#                      correction=c("border", "isotropic", "Ripley", "translate"),
#                     nlarge=Inf, domain=NULL, var.approx=FALSE, ratio=FALSE) {
  
  # Error checking
  verifyclass(X,"ppp")
  nslaves <- mpi.comm.size(comm) - 1
	if(verbose) print(paste(nslaves,"slaves found running on comm",comm))
  if (nslaves < 1) stop("Thre are no slaves running.")
  if (job.num < 2) stop("job.num needs to be at least 2.")
  if (job.num < nslaves) stop("job.num needs equal or larger than number of running slaves")

	# Send the comm to the slaves
	mpi.bcast.Robj2slave(comm)
	mpi.bcast.Robj2slave(verbose)
	
  if(verbose) print("No errors found, sending slaves into loops") 

	# Have all slaves go into a loop to await commands from master
  mpi.bcast.cmd(slave.parKest(), comm=comm)

  # Send number of jobs to all slaves (used for tag tracking)
  mpi.bcast(as.integer(job.num), type=1, comm=comm)
  
  # Send the ppp and arguments to each slave
  mpi.bcast.Robj(list(X=X, dot.arg=list(...)), rank=0, comm=comm)
 
  if(verbose) print("Objects sent to slaves, awaiting OK status")
		
	# Receive an OK status from slaves before continuting
	for (i in 1:nslaves) {
		tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
		srctag <- mpi.get.sourcetag()

		if(verbose) print(paste("--Slave",srctag[1],"ready to receive jobs"))
	}

  if(verbose) print("Global information sent to slaves, dividing subwindows")
  
  # Create a list of all subwindows
  subwin <- as.list(integer(job.num)) #preallocate list to store subwindows
  numpointper <- numeric(job.num) # vector to store number of points per subwindow
  for (i in 1:job.num) {
    xr <- X$window$xrange
    subxr <- c(xr[1]+(i-1)*diff(xr)/job.num,xr[1]+i*diff(xr)/job.num)
    subyr <- X$window$yrange
    slice <- owin(xrange=subxr,yrange=subyr)
    subwin[[i]] <- intersect.owin(slice,X$window)
    numpointper[i] <- X[subwin[[i]]]$n
  }
  if(verbose) print("Subwindows divided, sending first jobs to slaves")
  
  # Sort the jobs by size in each window from highest to lowest
  job.order <- order(numpointper,decreasing=TRUE)
  
  if(verbose) print("Subwindows sorted by decreasing difficulty")
 
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
  
  return(out)
}

slave.parKest <- function() {
  if (!exists("verbose")) verbose=FALSE
  if(verbose) print("Slave initialized, awaiting global information") 

  # get total number of jobs from master
  job.num <- mpi.bcast(integer(1), type=1, rank=0, comm=1)
  
  # get ppp and other arguments from master
  tmp <- mpi.bcast.Robj(rank=0, comm=1)
  X <- tmp$X
  dotarg <- tmp$dot.arg

	# send ok status to master
	mpi.send.Robj(obj=1, dest=0, tag=1, comm=1)
  
  if(verbose) print("Global information received, entering loop to receive jobs")

  #continually get jobs until there are no jobs left (indicated by tag)
  repeat {
    subwin <- mpi.recv.Robj(source=0, tag=mpi.any.tag(), comm=1)
    tag <- mpi.get.sourcetag()[2]
    print(tag)
    if (tag == 0)
      break
    
    K <- try(do.call("Kest", c(list(X=X, domain=subwin, nlarge=Inf, ratio=TRUE), dotarg)))
    attr(attr(K,"fmla"),".Environment") <- emptyenv() #remove unefficiently serialized environment
    
    mpi.send.Robj(obj=K, dest=0, tag=tag, comm=1)
  }
  if(verbose) print("Jobs completed, no more jobs, exiting loop")

}

#Modified version of pool.rat from spatstat that pools sub K estimates of the same pattern into a single one
pool.K <- function (...) 
{
  if (!is.list(...)) {
    argh <- list(...)
  } else {
    argh <- c(...)
  }
  
  n <- narg <- length(argh)
  if (narg == 0) 
    return(NULL)
  if (narg == 1) 
    return(argh[[1]])
  
  template <- vanilla.fv(argh[[1]])
  Y <- lapply(argh, attr, which = "numerator")
  X <- lapply(argh, attr, which = "denominator")
  templateX <- vanilla.fv(X[[1]])
  templateY <- vanilla.fv(Y[[1]])
  Add <- function(A, B) {
    force(A)
    force(B)
    eval.fv(A + B)
  }
  sumX <- Reduce(Add, X)
  sumY <- Reduce(Add, Y)
  attributes(sumX) <- attributes(templateX)
  attributes(sumY) <- attributes(templateY)
  Ratio <- eval.fv(sumY/sumX)
  attributes(Ratio) <-  attributes(template)
  #Ratio <- prefixfv(Ratio, tagprefix = "pool", descprefix = "pooled ", 
  #                  lablprefix = "")
  #result <- Reduce(bind.fv, list(Ratio, Variance, hiCI, loCI))
  return(Ratio)
}


parLest <- function(X, ..., job.num=2*(mpi.comm.size(comm)-1), comm=1, verbose=FALSE) {
  K <- parKest(X, ..., job.num, comm, verbose)
  L <- eval.fv(sqrt(K/pi))

  #need to add support for rip and ls variance estimation
  L <- rebadge.fv(L, new.ylab=quote(L(r)), new.fname="L", tags=names(K), new.labl=attr(K,"labl"))
  return(L)
}
