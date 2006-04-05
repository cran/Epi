\name{detrend}
\alias{detrend}
\title{ Projection of a model matrix on to the orthogonal
  complement of a trend. }
\description{
  The columns of the model matrix \code{M} is projected on the
  orthogonal complement to the matrix \code{(1,t)}. Orthogonality 
  is defined w.r.t. an inner product defined by the weights \code{weight}.
}
\usage{
  detrend( M, t, weight = rep(1, nrow(M)) )
}
\arguments{
  \item{M}{ A model matrix. }
  \item{t}{ The trend defining a subspace. A numerical vector of length
    \code{nrow(M)} }
  \item{weight}{ Weights defining the inner product of vectors \code{x}
    and \code{y} as \code{sum(x*w*y)}.
    A numerical vector of length \code{nrow(M)}, defaults to a vector of
    \code{1}s.}
}
\details{
  The functions is intended to be used in parametrization of
  age-period-cohort models.
}
\value{
  A full-rank matrix with columns orthogonal to \code{(1,t)}.
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://www.pubhealth.ku.dk/~bxc}, with help from Peter Dalgaard.
}
\seealso{ \code{\link{projection.ip}} }
\keyword{array}
