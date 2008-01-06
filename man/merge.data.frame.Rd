\name{merge.data.frame}
\alias{merge.data.frame}
\title{Merge data frame with a Lexis object}
\description{
  Merge two data frames, or a data frame with a \code{Lexis} object.
}
\usage{
\method{merge}{data.frame}(x, y, ...)
}
\arguments{
  \item{x, y}{data frames, or objects to be coerced into one}
  \item{...}{optional arguments for the merge method}
}
\details{
  This version of \code{merge.default} masks the one in the \code{base}.
  It ensures that, if either \code{x} or \code{y} is a \code{Lexis}
  object, then \code{merge.Lexis} is called.
}
\value{
  A merged \code{Lexis} object or data frame.
}
\author{Martyn Plummer}
\seealso{\code{\link{Lexis}}}
\keyword{ts}
