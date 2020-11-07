\name{detrend}
\alias{detrend}
\alias{decurve}
\title{ Projection of a model matrix on the orthogonal
  complement of a trend or curvature.}
\description{
  The columns of a model matrix \code{M} is projected on the
  orthogonal complement to the matrix \code{(1,t)},
  resp. \code{(1,t,t^2)}.
  
  Orthogonality is w.r.t. an inner product defined by the positive
  definite matrix matrix \code{diag(weight)}. Non-diagonal matrices
  defining the inner product is not supported.
}
\usage{
  detrend( M, t, weight = rep(1, nrow(M)) )
  decurve( M, t, weight = rep(1, nrow(M)) )
}
\arguments{
  \item{M}{A model matrix.}
  \item{t}{The trend defining a subspace. A numerical vector of length
    \code{nrow(M)}.}
  \item{weight}{ Weights defining the inner product of vectors \code{x}
    and \code{y} as \code{sum(x*w*y)}.
    A numerical vector of length \code{nrow(M)}, defaults to a vector of
    \code{1}s. Must be all non-negative.}
}
\details{
  The functions are intended to be used in construction of particular
  parametrizations of age-period-cohort models. 
}
\value{
  \code{detrend} returns full-rank matrix with columns orthogonal to
  \code{(1,t)};
  \code{decurve} returns full-rank matrix with columns orthogonal to
  \code{(1,t,t^2)}.
}
\author{
  Bendix Carstensen, Steno Diabetes Center Copenhagen,
  \url{http://bendixcarstensen.com}, with essential help from Peter Dalgaard.
}
\seealso{ \code{\link{projection.ip}} }
\keyword{array}
