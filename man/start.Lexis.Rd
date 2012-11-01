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
 entry(x, time.scale = NULL, by.id=FALSE)
  exit(x, time.scale = NULL, by.id=FALSE)
status(x, at="exit"        , by.id=FALSE)
   dur(x,                    by.id=FALSE)
}
\arguments{
  \item{x}{an object of class \code{Lexis}.}
  \item{time.scale}{a string or integer indicating the time scale. If
    omitted, all times scales are used.}
  \item{by.id}{Logical, if \code{TRUE}, only one record per unique value
  of \code{lex.id} is returned; either the first, the last or for
  \code{dur}, the sum of \code{lex.dur}. If \code{TRUE}, the returned
  object have the \code{lex.id} as (row)nmes attribute.}
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
