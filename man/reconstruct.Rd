\name{reconstruct}
\alias{reconstruct}
%\alias{parreconstruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
%Parallel and non-parallel point pattern reconstruction based on matching summary characteristics.
Point pattern reconstruction based on matching summary characteristics.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a point pattern and a set of summary characteristics, a new point pattern is reconstructed to best match the summary characteristics computed on the original point pattern. The reconstructed pattern can be the same size as the original or can be extended to a larger extent.
%This algorithm is implemented in parallel as well to speed up convergence.
%% if parallel reading is supported, can have slaves read invidually too
}
\usage{
reconstruct(x,n,fun=c("Kest"),correction=c("isotropic"),distance="absolute",
        win,eps=0.01,m=100,maxiter=5000,verbose=FALSE,balance=FALSE,
        conditional=FALSE)
%parreconstruct(x,n,fun=c("Kest"),correction=c("isotropic"),distance="absolute",win,
%       eps=0.01,m=100,maxiter=5000,verbose=FALSE,balance=FALSE,nslaves=2,nsteps=5,conditional=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The point process pattern (ppp) object used as a target for reconstruction.
}
  \item{n}{
Optional. The number of points to reconstruct, this is typically not given and is instead computed from the intensity of \code{x} and the given spatial window \code{win}.
  }
  \item{fun}{
The function on which to target reconstruction on. If you wish to compare on multiple estimates, a function can be written to incorporate these estimates to return one value or vector on which a comparison can be done. Defaults to Ripley's K-function.
}
  \item{correction}{
The type of edge correction used for \code{fun}, this is only required if \code{fun} requires it. Defaults to isotropic correction (assuming Ripley's K-function).
}
  \item{distance}{
How distance between two summary functions should be measured. Currently only absolute difference is supported.
}
  \item{win}{
The spatial window of class \code{\link{owin}} on which a point pattern will be reconstructed. If no spatial window is given, the spatial window for \code{x} will be used.
}
  \item{eps}{
The maximum improvement to be considered not having enough improvement.
}
  \item{m}{
The number of iterations without reaching \code{eps} amount of improvement before considered converge.
}
  \item{maxiter}{
The maximum number of iterations to run. The actual number of iterations may be lower if the stopping criteria is reached (\code{m} steps without at least \code{eps} improvement).
}
  \item{verbose}{
Optional. Whether or not each iteration with an improvement should be outputted to screen along with other debugging output.  
}
  \item{balance}{
Optional. Whether or not to use a balanced reconstruction algorithm to not favour long-range behavior of summary functions.
}
  \item{conditional}{
If reconstructing a pattern to a larger extent, this specifies whether or not the original pattern should be held fixed in the larger extent. No new points will be introduced to the original window.
}
%  \item{nslaves}{
%\code{parreconstruct} only. When reconstructing in parallel, the number of workers to use in parallel.
%}
%  \item{nsteps}{
%\code{parreconstruct} only. When reconstructing in parallel, the number of steps each worker should take at a time. The more %steps, the less communication but improvement will also diminish for higher number of steps.
%}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{ppp}{
An object of class \code{\link{ppp}} that contains the reconstructed point pattern.
}
\item{energy}{
A vector of the energy at each iteration.
}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Jonathan Lee \email{jlee253@uwo.ca}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
  # Reconstruct a cluster pattern based on K-function
  X <- rMatClust(12, 0.04, 7)
  Y <- reconstruct(X,eps=0.0001)
  par(mfrow=c(1,2))
  plot(X)
  plot(Y$ppp)
%  \dontrun{# Repeat the same thing but in parallel}
%  \dontrun{mpi.spawn.Rslaves(nslaves=4)}
%  \dontrun{mpi.bcast.cmd(require(parspatstat))}
%  \dontrun{X <- rMatClust(12, 0.04, 7)}
%  \dontrun{Y <- parreconstruct(X, eps=0.0001, nslaves=4, nsteps=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
