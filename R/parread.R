# read in a file and distribute to all slaves
# file = the filename to read from
# win = the window, if none then defaults to bounding box
# header = whether or not header row is present
# sep = separate input (defaults to comma)
# chunksize = the number of rows to read at a time
# xy = the column number of x and y coordinates
# marks = the column number of the marks for (x,y)
# xslice = whether the window should be sliced horizontally or vertically
# ... other arguments to pass to read.table
# job.num = number of jobs
# comm = communicator to use
# verbose = debug output

# reads it into a ppp object called localname locally on slaves

# returns a parppp object describing how the ppp is split amongst slaves, can be read by other parspatstat functions

parread.ppp <- function(file, win, header=TRUE, sep=",", chunksize=1000, xy=c(1,2), marks=NULL, xslice=TRUE, ..., job.num=(mpi.comm.size()-1), comm=1, verbose=FALSE, localname="X.ppp") {
  
  # Error checking
  if (missing(file)) {
    stop("Error: no file supplied. See help for examples.")
  }
    
  nslaves <- mpi.comm.size(comm) - 1
  if(verbose) print(paste(nslaves,"slaves found running on comm",comm))
  if (nslaves < 1) stop("Error: There are no slaves running.")
  if (job.num < 2) stop("Error: job.num needs to be at least 2.")
  if (job.num != nslaves) stop("Error: job.num needs to be equal to the number of running slaves")
  
  
  # if no window is defined, then go through data to find bounding window
  if (missing(win)) {
    if(verbose) print(paste("-No spatial window defined, finding bounding box"))
    xrange <- c(NA,NA)
    yrange <- c(NA,NA)
    eof <- FALSE
    numchunks <- 0
    while (!eof) {
      # read in chunksize+1 to account for boundary case where rows%%chunksize=0
      chunk <- read.table(file, nrows=chunksize+1, skip=(chunksize*numchunks)+header, sep=sep, header=FALSE)[,c(xy,marks)]
      
      numchunks <- numchunks + 1 # increment number of chunks read

      xrange <- c(min(xrange,min(chunk[-(chunksize+1),1]),na.rm=TRUE),max(xrange,max(chunk[-(chunksize+1),1]),na.rm=TRUE))
      yrange <- c(min(yrange,min(chunk[-(chunksize+1),2]),na.rm=TRUE),max(yrange,max(chunk[-(chunksize+1),2]),na.rm=TRUE))
            
      if (dim(chunk)[1] < (chunksize+1)) eof <- TRUE
    }
    
    win <- bounding.box.xy(xrange,yrange)
    if(verbose) print(paste("-Using bounding box as spatial window"))
    print(win)
  }
  
  #divide the window into subwindows
  subwin <- owin.divide(win,job.num,xslice=xslice)
  if(verbose) print(paste("-Spatial window divided into ",job.num, " strips."))
  
  # Have all slaves go into a loop to await commands from master
  #mpi.bcast.cmd(X <- slave.parread.ppp(), comm=comm)
  #mpi.remote.exec(cmd=slave.parread.ppp,verbose=verbose,localname=localname,comm=comm)
  mpi.bcast.cmd(cmd=assign,x=localname,value=0)
  mpi.bcast.cmd(cmd=assign,x=paste(localname,".subwin",sep=""),value=0)
  mpi.bcast.cmd(cmd=slave.parread.ppp,verbose=verbose,localname=localname,communicator=comm,comm=comm)
  # can't use remote.exec because it'll wait for a return before continuing
  
  if(verbose) print(paste("-Sending all slaves into loop to receive chunks."))

  # Send spatial window information to all slaves
  mpi.bcast.Robj(win, rank=0, comm=comm) #send full window
  for (i in 1:nslaves) {
    mpi.send.Robj(subwin[[i]], dest=i, tag=1, comm=comm) # send subwindow
  }
  if(verbose) print(paste("-Spatial window information sent to all slaves."))
  
  # Receive an OK status from slaves before continuting
  for (i in 1:nslaves) {
    tmp <- mpi.recv.Robj(source=mpi.any.source(), tag=mpi.any.tag(), comm=comm)
    srctag <- mpi.get.sourcetag()
    
    if(verbose) print(paste("--Slave",srctag[1],"ready to receive chunks"))
  }
  
  #only send points to the slaves that need it
  numchunks <- 0
  eof <- FALSE
  while(!eof) {  
    # read in chunksize+1 to account for boundary case where rows%%chunksize=0
    chunk <- read.table(file, nrows=chunksize+1, skip=(chunksize*numchunks)+header, sep=sep, header=FALSE)[,c(xy,marks)]
    
    numchunks <- numchunks + 1 # increment number of chunks read
    
    # check to see which slave the point should go to
    for (i in 1:nslaves) {
      ind <- inside.owin(chunk[-(chunksize+1),1],chunk[-(chunksize+1),2],subwin[[i]])
      
      #make sure index not empty
      if (sum(ind) > 0) {
        #send the rows that lie within subwindow i, ignoring the extra row
        mpi.send.Robj(chunk[-(chunksize+1),][ind,],dest=i,tag=1,comm=comm)
      }
    }
          
    if (dim(chunk)[1] < (chunksize+1)) eof <- TRUE
  }
  if(verbose) print(paste("-End of file reached, read in ", numchunks, " chunks."))
  
  #send stop messages to all slaves
  for (i in 1:nslaves) {
    mpi.send.Robj(obj=1,dest=i,tag=0,comm=comm)
  }
  
  if(verbose) print(paste("Stop signal sent to all slaves, slaves should have their respective local ppp objects created as X."))
  
  # Return object of class par.ppp
  result <- new("parppp",localname=localname,nslaves=nslaves,comm=comm)
  return(result)
}

# function that slaves run to receive data points from parread.ppp
slave.parread.ppp <- function(verbose=TRUE,localname,communicator) {
  #if (!exists("verbose")) verbose=TRUE
  if(verbose) print("Slave initialized, awaiting global information") 
  
  # receive full spatial window from master
  win <- mpi.bcast.Robj(rank=0, comm=communicator)
  if(verbose) print(win)
  
  # receive subwindow from master
  subwin <- mpi.recv.Robj(source=0, tag=mpi.any.tag(), comm=communicator)
  
  # send ok status to master
  mpi.send.Robj(obj=1, dest=0, tag=1, comm=communicator)
  
  if(verbose) print("Global information received, entering loop to receive jobs")
  
  #continually get jobs until there are no jobs left (indicated by tag)
  chunks <- list()
  chunknum <- 1
  repeat {
    packet <- mpi.recv.Robj(source=0, tag=mpi.any.tag(), comm=communicator)
    tag <- mpi.get.sourcetag()[2]
    if (tag == 0) break
    
    chunks[[chunknum]] <- packet
     
    if(verbose) print(paste("Chunk received. Chunk number ",chunknum))
    chunknum <- chunknum+1
  }
  
  if(verbose) print("All chunks received, creating final ppp")
    
  x <- numeric(0)
  y <- numeric(0)
  z <- numeric(0)
  for (i in 1:length(chunks)) {
    x <- c(x,chunks[[i]][,1])
    y <- c(y,chunks[[i]][,2])
    if (dim(chunks[[1]])[2] > 2)  z <- c(z,chunks[[i]][,3])
  }
  
  marks = NULL
  if (length(z) > 0) marks = z
  
  final.ppp <- ppp(x=x,y=y,window=win,marks=marks)
  if(verbose) print(final.ppp)
  
  #objects are saved in base namespace since global will cause R 3.0.0+ to give a NOTE with --as-cran CHECK
  assign(x=localname,value=final.ppp,envir=.BaseNamespaceEnv)
  assign(x=paste(localname,".subwin",sep=""),value=subwin,envir=.BaseNamespaceEnv)
  
}






