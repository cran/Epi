\name{start.Lexis}
\alias{entry}
\alias{exit}
\alias{status}
\alias{dur}
\title{Time series methods for Lexis objects}
\description{
  Extract the entry time, exit time, status, or duration of follow-up from a
  \code{Lexis} object.
}
\usage{
entry(x, time.scale = NULL)
exit(x, time.scale = NULL)
status(x, at="exit")
dur(x)
}
\arguments{
  \item{x}{an object of class \code{Lexis}.}
  \item{time.scale}{a string or integer indicating the time scale. If
    omitted, all times scales are used.}
  \item{at}{string indicating the time point(s) at which status is to be
    measured.}
}
\value{
  The \code{entry} and \code{exit} functions return a vector of 
  entry times and exit times, respectively, on the requested time
  scale.  If multiple time scales are requested, then a matrix is 
  returned.

  The \code{status} function returns a vector giving the status
  at entry or exit and \code{dur} returns a vector with the lengths
  of the follow-up intervals.
}
\author{Martyn Plummer}
\seealso{\code{\link{Lexis}}}
\keyword{survival}
\keyword{ts}
