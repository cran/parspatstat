\name{parKest}
\alias{parKest}
\alias{parLest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Parallel K-function and L-function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimate's Ripley's reduced second moment functino (K-function) and the transformed L-function in parallel. Slaves assumed to be already spawned via \code{\link{mpi.spawn.Rslaves}}.
}
\usage{
parKest(X, ..., job.num = 2 * (mpi.comm.size(comm) - 1), comm = 1, verbose = FALSE)
parLest(X, ..., job.num = 2 * (mpi.comm.size(comm) - 1), comm = 1, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
A point pattern of class \code{\link{ppp}}, for which K(r) is to be computed.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
Extra arguments to be passed into the \code{\link{Kest}} function.
}
  \item{job.num}{
%%     ~~Describe \code{job.num} here~~
The number of subwindows to split the point pattern into. This must be at least the number of slaves spawned and recommended (and defaults) to be at least twice as much for load balancing. If you are finding that each job is still too large to run on a single machine, try increasing this value.
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
The communicator the slaves are on. Defaults to 1.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Optional. Whether or not to print debugging output. Defaults to FALSE.
}
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
See \code{\link{Kest}} for details of value returned.
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{Kest}}
}
\examples{
  # Generate a random point pattern of 10000 points and compute it's K-estimate using 4 slaves in parallel.
  \dontrun{mpi.spawn.Rslaves(nslaves=4)}
  \dontrun{mpi.bcast.cmd(require(parspatstat))}
  \dontrun{X <- runifpoint(10000)}
  \dontrun{K <- parKest(X)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~spatial }
