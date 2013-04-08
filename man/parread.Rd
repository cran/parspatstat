\name{parread.ppp}
\alias{parread.ppp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Parallel reading of a file into a parallel point process pattern object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reads in a point process pattern and distributes appropriate pieces to slaves. Slaves assumed to be already spawned via \code{\link{mpi.spawn.Rslaves}}. Unique point patterns and subwindows will exist within each slave. Details of what is stored on each slave is returned by this function to the master.
%% if parallel reading is supported, can have slaves read invidually too
}
\usage{
parread.ppp(file, win, header=TRUE, sep=",", chunksize=1000, xy=c(1,2), marks=NULL, xslice=TRUE, ..., job.num=(mpi.comm.size()-1), comm=1, verbose=FALSE, localname="X.ppp")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The file from which a point pattern is to be read from.
}
  \item{win}{
The spatial window of class \code{\link{owin}} on which the point pattern is to exist. If no spatial window is provided, the smallest bounding box of all the point data will be used.
}
  \item{header}{
Whether or not the first row of the file contains names of the columns. Defaults to \code{TRUE}.
}
  \item{sep}{
The delimiter of the file, defaults to "," for comma separated values (.csv). Other common choices are space delimited or tab delimited.
}
  \item{chunksize}{
The number of rows for the master to read at a time. This value defaults to 1000 but can be increased to approach the amount of memory available to the master.  A smaller value is fine but will require more passes through the file to read it completely.
}
  \item{xy}{
A vector of length two indicating the two column numbers to read the x- and y-coordinate from respectively. Defaults to the first and second column.
}
  \item{marks}{
Optional. Indicates the column number of the marks (to be attached to the coordinate of the same row), defaults to NULL (unmarked point pattern).
}
  \item{xslice}{
Whether or not data should be partitioned horizontally or vertically. Defaults to TRUE (horizontal partitions).
}
  \item{\dots}{
Additional arguments to be passed into \code{read.table}.  
}
  \item{job.num}{
%%     ~~Describe \code{job.num} here~~
The number of subwindows to split the point pattern into. This must be at least the number of slaves spawned and recommended (and defaults) to be at least twice as much for load balancing. If you are finding that each job is still too large to run on a single machine, try increasing this value.
Important note: If the data is not pre-sorted, \code{job.num} must be exactly equal to the number of slaves spawned since load balancing will not be supported.
}
  \item{comm}{
%%     ~~Describe \code{comm} here~~
The communicator the slaves are on. Defaults to 1.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Optional. Whether or not to print debugging output. Defaults to FALSE.
}
  \item{localname}{
%%     ~~Describe \code{localname} here~~
The variable name under which each slave should save its ppp object. Subwindows are also saved on each slave with the format \code{localname.subwin}. Defaults to "X.ppp" but should be changed if more than one \code{\link{parppp}} object is to be saved on each slave.
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
An object of class \code{\link{parppp}} is returned. Each slave will contain an object named \code{localname} and \code{localname.subwin} that contains the point pattern and spatial window it is responsible for.
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
  \dontrun{mpi.spawn.Rslaves(nslaves=4)}
  \dontrun{mpi.bcast.cmd(require(parspatstat))}
  \dontrun{P <- parread.ppp(file="lightning.csv",chunksize=1000, xy=c(5,6), header=TRUE, sep=",", xslice=TRUE)}
  \dontrun{K.all <- parKest(X=P, verbose=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
