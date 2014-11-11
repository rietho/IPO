\name{sendRetGroupSlaveFunctions}
\alias{sendRetGroupSlaveFunctions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sendRetGroupSlaveFunctions(parameters, xset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (parameters, xset) 
{
    mpi.bcast.Robj2slave(calculateRGTV)
    mpi.bcast.cmd(slave <- mpi.comm.rank())
    mpi.bcast.Robj2slave(parameters)
    mpi.bcast.Robj2slave(xset)
    mpi.bcast.cmd("library(xcms)")
    mpi.bcast.Robj2slave(optimizeRetGroupSlave)
    mpi.bcast.cmd(optimizeRetGroupSlave())
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
