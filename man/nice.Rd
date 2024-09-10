\name{nice}
\alias{nice}
\title{Nice breakpoints for axes on plots}
\description{The function calls \code{\link{pretty}}
  for linear scale. For a log-scale nice are computed using a set of
  specified number in each decade.
}
\usage{
nice(x, log = FALSE, lpos = c(1, 2, 5), xmx = 4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numerical vector to}
  \item{log}{Logical. Is the scale logartimic?}
  \item{lpos}{Numeric. Numbers between 1 and 10 giving the desired
    breakpoints in this interval.}
  \item{xmx}{Numeric. The maximal (absolute) power of 10 to be used for
    a log-scale.}
  \item{\dots}{Arguments passed on to \code{pretty} if \code{log}=FALSE}
}
\value{A vector of breakpoints.}
\author{Bendix Carstensen, \email{b@bxc.dk},
  \url{http://bendixcarstensen.com}}
\seealso{pretty}
\examples{
nice( exp( rnorm( 100 ) ), log=TRUE )
}
\keyword{manip}
