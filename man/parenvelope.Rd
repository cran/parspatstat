\name{parenvelope}
\alias{parenvelope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Parallel Simulation Envelopes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes simulation envelopes of a summary function in parallel. Slaves assumed to be already spawned via \code{\link{mpi.spawn.Rslaves}}. Usage is almost identical to the \code{\link{envelope}} function.
}
\usage{
parenvelope(Y, fun=Kest, nsim=99, ..., job.num=2*(mpi.comm.size(comm)-1), 
                 verbose=FALSE, comm=1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Either a point pattern (\code{\link{ppp}}) object or a fitted point process model (\code{\link{ppm}}) on which the summary statistic will be computed.
}
  \item{fun}{
%%     ~~Describe \code{\dots} here~~
The summary statistic function on which the envelope is simulated. Parallel summary functions (e.g. \code{\link{parKest}}, \code{\link{parLest}}) can be used as well, but slaves need to be spawned on a separate communicator and specified accordingly.
}
  \item{nsim}{
%%     ~~Describe \code{\dots} here~~
Number of total simulations to be run -- will be split amongst slaves.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
Extra arguments to be passed into the \code{\link{envelope}} function.
}
  \item{job.num}{
%%     ~~Describe \code{comm} here~~
The number of jobs to split the simulations into for load balancing. Defaults to two times the number of slaves.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Optional. Whether or not to print debugging output. Defaults to FALSE.
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
The communicator the slaves are on. Defaults to 1.
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
\code{\link{envelope}}
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
\code{\link{envelope}}
}
\examples{
  # Generate a clustered point pattern of ~4000 points and generate a simulation 
  # envelope of the K-function with 1000 simulations in parallel.
  \dontrun{mpi.spawn.Rslaves(nslaves=4)}
  \dontrun{mpi.bcast.cmd(require(parspatstat))}
  \dontrun{X <- rMatClust(1000, 0.05, 4)}
  \dontrun{K <- parenvelope(X, fun=Kest, nsim=1000)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
