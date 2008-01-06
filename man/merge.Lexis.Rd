\name{merge.Lexis}
\alias{merge.Lexis}
\title{Merge a Lexis object with a data frame}
\description{
  Merge additional variables from a data frame into a Lexis object.
}
\usage{
\method{merge}{Lexis}(x, y, id, by, ...)
}
\arguments{
  \item{x}{an object of class \code{Lexis}}
  \item{y}{a data frame}
  \item{id}{the name of the variable in \code{y} to use for matching
    against the variable \code{lex.id} in \code{x}.
  }
  \item{by}{if matching is not done by id, a vector of variable names
    common to both \code{x} and \code{y}}
  \item{...}{optional arguments to be passed to \code{merge.data.frame}}
}
\details{
  A \code{Lexis} object can be considered as an augmented data frame
  in which some variables are time-dependent variables representing
  follow-up. The \code{Lexis} function produces a minimal object
  containing only these time-dependent variables.  Additional variables
  may be added to a \code{Lexis} object using the \code{merge} method.
}
\value{
  A \code{Lexis} object with additional columns taken from the
  merged data frame.
}
\author{Martyn Plummer}
\note{
  The variable given as the \code{by.y} argument must not contain
  any duplicate values in the data frame \code{y}.
  }
\seealso{\code{\link{merge.data.frame}}, \code{\link{subset.Lexis}}}
\keyword{array}
\keyword{manip}
