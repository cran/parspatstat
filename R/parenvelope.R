parenvelope <- function(Y, fun=Kest, nsim=99, ..., job.num=2*(mpi.comm.size(comm)-1), verbose=FALSE, comm=1) {
  
  # Error checking
  nslaves <- mpi.comm.size(comm) - 1
  if(verbose) print(paste(nslaves,"slaves found running on comm",comm))
  if (nslaves < 1) stop("There are no slaves running.")
  if (job.num < 2) stop("job.num needs to be at least 2.")
  if (job.num < nslaves) stop("job.num needs equal or larger than number of running slaves")
  if (nsim < nslaves) stop("Cannot run less simulations than there are slaves.")
  
  mpi.bcast.Robj2slave(verbose)
  
  if(verbose) print("No errors found, sending slaves into function") 
  
  # Have all slaves go into a function to receive the following commands from master
  mpi.bcast.cmd(cmd=slave.parenvelope, verbose=verbose, job.num=job.num, Y=Y, fun=fun, ..., communicator=comm, comm=comm)

  # Send number of jobs to all slaves (used for tag tracking)
  #mpi.bcast(as.integer(job.num), type=1, comm=comm)
  
  # Send Y and arguments to each slave
  #mpi.bcast.Robj(list(Y=Y, fun=fun, dot.arg=list(...)), rank=0, comm=comm)
  
  if(verbose) print("Objects sent to slaves, awaiting OK status")
  
  # Receive an OK status from slaves before continuting
  for (i in 1:nslaves) {
    tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
    srctag <- mpi.get.sourcetag()
    
    if(verbose) print(paste("--Slave",srctag[1],"ready to receive job"))
  }
  
  # Split nsim into the number of jobs we have
  jobs <- rep(floor(nsim/job.num),job.num)
  if (nsim%%job.num > 0) jobs[1:(nsim%%job.num)] <- jobs[1]+1
  
  if(verbose) print("nsim divided, sending each nsim to slaves")
  
  # Send nsim to slaves
  #for (i in 1:nslaves) {
  #  mpi.send.Robj(jobs[i], dest=i, tag=i, comm=comm)
  #  if(verbose) print(paste("--Job",i,"sent to slave", i))
  #}
  
  # Send first nslaves job to slaves
  for (i in 1:nslaves) {
    mpi.send.Robj(jobs[i], dest=i, tag=i, comm=comm)
    if(verbose) print(paste("--Job",i,"sent to slave", i))
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
      mpi.send.Robj(obj=jobs[j], dest=srctag[1], tag=j, comm=comm)
      if(verbose) print(paste("--Job",j,"sent to slave", srctag[1]))
    } else {
      # otherwise tell the slave it's finished
      mpi.send.Robj(obj=as.integer(0), dest=srctag[1], tag=0, comm=comm)
      if(verbose) print(paste("--Termination command sent to slave", srctag[1]))
    }
  }
  
  if(verbose) print("All jobs complete and received, formatting output")
  
  # Combine the simulated envelope objects into a single one and return it
  Eobj <- do.call(what=pool.envelope,args=out)
  
  return(Eobj)
  
}


slave.parenvelope <- function(verbose=FALSE, job.num, Y, fun, ..., communicator) {
  #if (!exists("verbose")) verbose=FALSE
  if(verbose) print("Slave initialized, awaiting global information") 
  
  # get total number of jobs from master
  #job.num <- mpi.bcast(integer(1), type=1, rank=0, comm=1)
  
  # get Y and other arguments from master
  #tmp <- mpi.bcast.Robj(rank=0, comm=1)
  #Y <- tmp$Y
  #fun <- tmp$fun
  #dotarg <- tmp$dot.arg
  
  # send ok status to master
  mpi.send.Robj(obj=1, dest=0, tag=1, comm=communicator)
  
  if(verbose) print("Global information received, retrieving number of simulations")
  
  #continually get jobs until there are no jobs left (indicated by tag)
  repeat {
    nsim <- mpi.recv.Robj(source=0, tag=mpi.any.tag(), comm=communicator)
    tag <- mpi.get.sourcetag()[2]
    print(tag)
    if (tag == 0)
      break
    
    output <- envelope(Y=Y,fun=fun,nsim=nsim,savefuns=TRUE,...)
    #output <- try(do.call(envelope, c(list(Y=Y, fun=fun, nsim=nsim, savefuns=TRUE), dotarg)))

    mpi.send.Robj(obj=output, dest=0, tag=tag, comm=communicator)
  }
    
  if(verbose) print("Simulation complete")
}

