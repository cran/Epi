\name{subset.Lexis}
\alias{subset.Lexis}
\title{Subsetting Lexis objects}
\description{
  Return subsets of Lexis objects which meet conditions
}
\usage{
\method{subset}{Lexis}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{Lexis}}
  \item{\dots}{additional arguments to be passed to \code{subset.data.frame}}
}
\details{
  The subset method for \code{Lexis} objects works exactly as the
  method for data frames.
}
\value{
  A \code{Lexis} object with selected rows and columns.
}
\author{Martyn Plummer}
\seealso{\code{\link{Lexis}}, \code{\link{merge.Lexis}}}
\keyword{manip}
