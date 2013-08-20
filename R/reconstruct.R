#reconstruction algorithm based on a given statistic (or combination of statistics)

# x is the target point pattern (should be ppp object)
# fun is a list of target summary functions, if more than one, sum is taken
#   Kest, Gest, Hest, Jest, Lest
# correction is a list of correction methods to use from the specified summary function
#   in the same order as 'target'
# ... additional arguments passed into target functions
# distance is distance measure, defaults to absolute distance
# win is the target spatial window on which to reconstruct the point pattern
# error is tolerance before stopping
reconstruct <- function(x,n,fun=c("Kest"),correction=c("isotropic"),distance="absolute",win,eps=0.01,m=100,maxiter=5000,verbose=FALSE,balance=FALSE,conditional=FALSE) {
    
  targetfuns <- list()
  # compute target summaries and save into a list
  for (i in 1:length(fun)) {
    Fx <- do.call(fun[i],args=list(X=x,correction=correction[i]))
    targetfuns <- c(targetfuns,list(Fx))
  }
  
  init.win = x$win
  if (!missing(win)) {
    init.win = win
  }
  
  init.n = x$n
  if (!missing(n)) {
    init.n = n
  } else if (missing(n) & missing(win)) {
    init.n = x$n
  } else if (missing(n) & !missing(win)) {
    init.n = ceiling(n/area.owin(x$win)*area.owin(win)) #scale n up (or down) to new window size
  }
  # generate initial random pattern
  y <- runifpoint(init.n,win=init.win)
  
  # fix the point pattern to include conditional points
  if (conditional) {
    tempy <- runifpoint(init.n+500,win=init.win) #some extra points
    tempy <- tempy[-ppp.extract(tempy,x$win)] #cut out small window
    
    # strip excess points
    if (tempy$n > init.n - x$n) {
      tempy <- tempy[1:(init.n-x$n)] 
    }
    
    y <- superimpose(tempy,x) #add in conditional window
  }
  
  # compute weight on distances based on radiuses
  weight = 1
  if (balance==TRUE) {
    weight = poissonVar(init.n,rep=500)+0.000001
  }
  
  
  # compute energy between current and target
  energy <- Inf
  iter <- 1
  while (iter < maxiter) {
    # randomly move a point in y
    newy <- y
    
    # set to only move points outside conditional area
    samplerange <- 1:newy$n
    if (conditional) {
      samplerange <- (1:newy$n)[-ppp.extract(newy,x$win)]
    }
    
    rand.i <- sample(samplerange,1)
    rand.x <- runif(1,min=newy$win$xrange[1],max=newy$win$xrange[2])
    rand.y <- runif(1,min=newy$win$yrange[1],max=newy$win$yrange[2])
    
    # resample point if it is inside the conditional window
    if (conditional) {
      while (length(ppp.extract(ppp(rand.x,rand.y,window=init.win),x$win))>0) {
        rand.x <- runif(1,min=newy$win$xrange[1],max=newy$win$xrange[2])
        rand.y <- runif(1,min=newy$win$yrange[1],max=newy$win$yrange[2])
      }
    }
    # check if point is actually in required spatial window
    # reselect otherwise
    
    newy$x[rand.i] <- rand.x
    newy$y[rand.i] <- rand.y
    
    newenergy <- 0
    for (i in 1:length(fun)) {
      Fy <- do.call(fun[i],args=list(X=newy,r=targetfuns[[i]]$r,correction=correction[i]))
      newenergy <- newenergy + sum((((targetfuns[[i]]$iso - Fy$iso)^2)/weight))
    } #end fun loop
    
    #check if newenergy is better
    if (newenergy < energy[iter]) {
      energy <- c(energy,newenergy)
      y <- newy
      if (verbose) print(paste("Iteration: ",iter,"; Energy: ",energy[iter],sep=""))
    } else {
      # keep same as last step
      energy <- c(energy,energy[iter])
    }
    
    # check if there has been less than eps improvement in last m iterations
    if (iter > m) {
      if (energy[iter-m] - energy[iter] < eps) {
        break
      }
    }
    #increment iteration count
    iter <- iter+1
  }
  
  if (verbose) {
    if (iter >= maxiter) print("Maximum number of iterations reached, reconstruction halted.")
    else print("Not enough improvement made in last m steps, halting.")
  }
  # return resulting point pattern
  return(list(ppp=y,energy=energy))
}




# run 10 sequentially
# discuss why krigingis not the same thing
# add indicator function for hard-core inhibition
# discuss non-blocking
# discuss windows with holes -- polygon vs mask


# normalize statistics


# default epsilon - something on scale of square root of variance
#   take largest variance
# change scale -- same model and # of points. stopping time averaged over many simulations should be about equal
# 


# compute K for smaller region and reconstruct a large region
# then compare with a larger test region

# can compare on other statistics

# model checking
# given a large dataset with Khat
# thin the large dataset to something smaller
# reconstruct the larger dataset with Khat to see if it matches the original



# fault tolerance -- PVM had it, MPI decided not to use it





