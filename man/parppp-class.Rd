\name{parppp-class}
\Rdversion{1.1}
\docType{class}
\alias{parppp-class}
\alias{print.parppp}
\alias{plot.parppp}
\alias{parppp}

\title{Class \code{"parppp"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
A class \code{parppp} to define a parallel point process pattern object that exists locally on each slave
}
\section{Objects from the Class}{
In the \pkg{parspatstat} library, a large point pattern can be read by a master and split amongst the slaves, so that each slave only contains a managable portion of the entire point pattern.  Each slave is still aware of the extent of the full window as well as the extent of the subwindow it is responsible for. 

  Objects of class \code{"parppp"} are created by the function
  \code{\link{parread.ppp}}

  Functions that can accept \code{"parppp"} objects include
  \code{\link{parKest}}
  and
  \code{\link{parLest}}.

  There is a \code{plot} method for parallel point pattern objects,
  \code{\link{plot.parppp}}. There is also a method for \code{print}. 
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{localname}:}{ the name under which a point pattern is saved locally on slaves.}
    \item{\code{nslaves}:}{the number of slaves that are involved in storing this parallel point pattern.}
    \item{\code{comm}:}{the communicator on which the slaves are on.}
  }
}
\section{Methods}{
No methods defined with class "parppp" in the signature.
}
\author{
%%  ~~who you are~~
Jonathan Lee \email{jlee253@uwo.ca}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{parread.ppp}},
  \code{\link{parKest}},
  \code{\link{parLest}}
}
\examples{
showClass("parppp")
}
\keyword{classes,spatial}
