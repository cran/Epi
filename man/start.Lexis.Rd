\name{start.Lexis}
\alias{entry}
\alias{exit}
\alias{status}
\alias{dur}
\alias{transient}
\alias{absorbing}
\alias{before}
\alias{preceding}
\alias{after}
\alias{succeeding}
\title{Time series and other methods for Lexis objects}
\description{
  Extract the entry time, exit time, status or duration of follow-up from a
  \code{Lexis} object. Classify states.
}
\usage{
     entry( x, time.scale = NULL, by.id=FALSE )
      exit( x, time.scale = NULL, by.id=FALSE )
    status( x, at="exit"        , by.id=FALSE )
       dur( x,                    by.id=FALSE )
 transient( x )
 absorbing( x )
 preceding( x, states )
    before( x, states )
succeeding( x, states )
     after( x, states )
}
\arguments{
  \item{x}{an object of class \code{Lexis}.}
  \item{time.scale}{a string or integer indicating the time scale. If
    omitted, all times scales are used.}
  \item{by.id}{Logical, if \code{TRUE}, only one record per unique value
  of \code{lex.id} is returned; either the first, the last, or for
  \code{dur}, the sum of \code{lex.dur}. If \code{TRUE}, the returned
  object have the \code{lex.id} as (row)names attribute.}
  \item{at}{string indicating the time point(s) at which status is to be
    measured. Possible values are "exit" or "entry".}
  \item{states}{Character vector of states.}
}
\value{
  The \code{entry} and \code{exit} functions return a vector of
  entry times and exit times, respectively, on the requested time
  scale.  If multiple time scales are requested, a matrix is
  returned.

  The \code{status} function returns a vector giving the status at
  "\code{at}" (entry or exit) and \code{dur} returns a vector with the
  lengths of the follow-up intervals.

  The functions \code{transient} and \code{absorbing} return character
  vectors of the transient, resp. absorbing states. These are necessariy
  disjoint but the union may not be the same as \code{levels(x)}, since
  the latter may have levels that are never assumed by either
  \code{lex.Cst} or \code{lex.Xst}.

  \code{preceding} returns a character vector with names of the states of the
  Lexis object \code{x} from which one of the states in
  \code{states} can be reached directly - the preceding
  states. \code{before} is just a synonym for \code{preceding}.

  \code{succeeding} returns a character vector with names of the states of the
  Lexis object \code{x} that can be reached directly from one of the
  states in \code{states}. \code{after} is just a synonym for \code{succeeding}.
}
\author{Martyn Plummer & Bendix Carstensen}
\seealso{\code{\link{Lexis}}}
\keyword{survival}
\keyword{ts}
