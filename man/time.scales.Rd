\name{timeScales}
\alias{timeScales}
\title{The time scales of a Lexis object}
\description{
  Function to get the names of the time scales of a \code{Lexis} object.
}
\usage{
timeScales(x)
}
\arguments{
  \item{x}{an object of class \code{Lexis}.}
}
\value{
  A character vector containing the names of the variables in \code{x}
  that represent the time scales. Extracted from the \code{time.scales}
  attribute of the object.
}
\author{Martyn Plummer}
\seealso{\code{\link{Lexis}}, \code{\link{splitLexis}}}
\keyword{attribute}
