\name{timeScales}
\alias{timeScales}
\alias{timeSince}
\alias{tsNA20}
\title{The time scales of a Lexis object}
\description{
  Functions to get the names and type of the time scales of a
  \code{Lexis} object. 
}
\usage{
timeScales(x)
timeSince(x)
tsNA20( x, all.scales=FALSE )
}
\arguments{
  \item{x}{an object of class \code{Lexis}.}
  \item{all.scales}{Should NAs in all timescales be replaced by 0? If
    \code{FALSE} (the default) only timescales defined as time since
    entry to a state get \code{NA}s replaced by 0s}
}
\value{
  \code{timeScales} returns a character vector containing the names of
  the variables in \code{x} that represent the time scales. Extracted
  from the \code{time.scales} attribute of the object.
  
  \code{timeSince} returns a named character vector, the names being the
  names of the timescales and the content being the names of the states
  to which the corresponding timescale is defined as time since
  entry. For those time scales that are not defined as such an empty
  string is used. Hence, if none of the timescales are defined as time
  since entry to a state \code{timeSince} will return a vector of empty
  strings.
}
\author{Martyn Plummer, Bendix Carstensen}
\seealso{\code{\link{Lexis}}, \code{\link{splitLexis}}}
\keyword{attribute}
