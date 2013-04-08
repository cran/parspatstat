parLest <- function(X, ...) {
  K <- parKest(X=X, ...)
  L <- eval.fv(sqrt(K/pi))
  
  #need to add support for rip and ls variance estimation
  L <- rebadge.fv(L, new.ylab=quote(L(r)), new.fname="L", tags=names(K), new.labl=attr(K,"labl"))
  return(L)
}