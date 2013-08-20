#parspatstat.R
# various helper functions for parspatstat

#function to sample uniformly from a point process pattern object
#different from rthin in that you specify the exact size
# returns ppp object
sample.ppp <- function(x,size,replace=FALSE,prob=NULL) {
  ind <- sample(1:x$n,size=size,replace=replace,prob=prob)
  return(x[ind])
}


# owin.divide - divide a spatial window into n pieces and return a list of subwindows
# win: the spatial window (owin object)
# n: number of slices
# xslice: slice in x-axis or y-axis
# R: (optional), expand the window by 0 units in one direction (depending on xslice)
owin.divide <- function(w, n, xslice=TRUE, R=0) {
  subwin <- as.list(integer(n)) #preallocate list to store subwindows
  
  if (xslice) { #slice in x axis
    for (i in 1:n) {
      xr <- w$xrange
      subxr <- c(xr[1]+(i-1)*diff(xr)/n,xr[1]+i*diff(xr)/n)
      subyr <- w$yrange
      slice <- owin(xrange=subxr,yrange=subyr)
      subwin[[i]] <- intersect.owin(slice,w)
    }
  } else {
    for (i in 1:n) {
      yr <- w$yrange
      subyr <- c(yr[1]+(i-1)*diff(yr)/n,yr[1]+i*diff(yr)/n)
      subxr <- w$xrange
      slice <- owin(xrange=subxr,yrange=subyr)
      subwin[[i]] <- intersect.owin(slice,w)
    }
  }
  
  return(subwin)
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

poissonVar <- function(n,rep=50) {
  Ks <- numeric(0)
  for (i in 1:rep) {
    X <- runifpoint(n)
    K <- Kest(X,correction="isotropic")
    Ks <- cbind(Ks,K$iso)
  }
  
  Ks.var <- apply(Ks,1,var)
  return(Ks.var)
}

# function that returns indices of x that are within owin 'win'
ppp.extract <- function(x,win) {
  ind <- which(x$x > win$xrange[1] & x$x < win$xrange[2] & x$y > win$yrange[1] & x$y < win$yrange[2])
  return(ind)  
}