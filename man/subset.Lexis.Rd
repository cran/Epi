\name{subset.Lexis}
\alias{subset.Lexis}
\alias{[.Lexis}
\alias{subset.stacked.Lexis}
\title{Subsetting Lexis (and stacked.Lexis) objects}
\description{
  Return subsets of Lexis objects which meet conditions
}
\usage{
\method{subset}{Lexis}(x, ...)
\method{[}{Lexis}(x, ...)
\method{subset}{stacked.Lexis}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{Lexis}}
  \item{\dots}{additional arguments to be passed to
    \code{subset.data.frame}. This will normally be some logical
    expression selecting a subset of the rows. For details see
    \code{\link{subset.data.frame}}.}
  }
\details{
  The subset method for \code{Lexis} objects works exactly as the method
  for data frames. So does the "[" method. The special methods are needed in
  order to propagate the Lexis-specific attributes.

  The method for \code{stacked.Lexis} objects also shrinks the set of
  levels for \code{lex.Cst} and \code{lex.Xst} to those actually
  occurring in the resulting data frame.
  }
\value{
  A \code{Lexis} object with selected rows and columns.
}
\author{Martyn Plummer}
\seealso{\code{\link{Lexis}}, \code{\link{merge.Lexis}}}
\keyword{manip}
