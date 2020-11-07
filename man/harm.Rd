\name{harm}
\alias{harm}
\title{Create a basis of harmonic functions.
}
\description{Returns a matrix of harmonic functions usable for modeling
periodic effects
}
\usage{
harm(x, ord=1, per=1, verbose=FALSE )
}
\arguments{
  \item{x}{A numeric variable.
}
  \item{ord}{Integer, the order of the harmonic.
}
  \item{per}{Numeric, the length of the period on the \code{x} scale.
}
  \item{verbose}{Logical: shall I tell what I do with dates?
}
}
\details{
  Columns are constructed under the assumption that the periodic function
  has period \code{per} on the \code{x} scale. Thus, the first columns
  is defined as \code{sin(2*pi*x/per)}, \code{cos(2*pi*x/per)},
  \code{sin(4*pi*x/per)} etc.

  Since \code{sin} and \code{cos} are periodic functions there is no
  requirement that \code{x} be in any particular range.
}
\value{A matrix with \code{nrow(x)} rows and
  2*\code{deg} columns and columnnames
  \code{sin1}, \code{cos1},
  \code{sin2}, \code{cos2} etc.
}
\author{
  Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\examples{
x <- seq(-1,1,0.01)
head( harm(x,ord=2) )
matplot( x, harm(x,ord=2), type="l", lty=1, lwd=3 )
}
\keyword{models}
\keyword{regression}
