\name{projection.ip}
\alias{projection.ip}
\title{ Projection of columns of a matrix. }
\description{
  Projects the columns of the matrix \code{M} on the space spanned by the
  columns of the matrix \code{X}, with respect to the inner product
  defined by \code{weight}: \code{<x|y>=sum(x*w*y)}.
  }
\usage{
projection.ip(X, M, orth = FALSE, weight = rep(1, nrow(X)))
}
\arguments{
  \item{X}{ Matrix defining the space to project onto. }
  \item{M}{ Matrix of columns to be projected. Must have the same number
    of rows as \code{X}. }
  \item{orth}{ Should the projection be on the orthogonal complement to
    \code{span(X)}? }
  \item{weight}{ Weights defining the inner product. Numerical vector of
  length \code{nrow(X)}. }
}
\value{
  A matrix of full rank with columns in \code{span(X)}.
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://BendixCarstensen.com}, with help from Peter Dalgaard.
}
\seealso{ \code{\link{detrend}} }
\keyword{array}
