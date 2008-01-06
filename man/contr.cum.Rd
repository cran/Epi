\name{contr.cum}
\alias{contr.cum}
\alias{contr.2nd}
\alias{contr.diff}
\alias{contr.orth}
\title{
  Contrast matrices
  }
\description{
  Return a matrix of contrasts for factor coding.
  }
\usage{
  contr.cum(n)
  contr.diff(n)
  contr.2nd(n)
  contr.orth(n) 
  }
\arguments{
  \item{n}{A vector of levels for a factor, or the number of levels.}
  }
\details{
  These functions are used for creating contrast matrices for use in
  fitting regression models. The columns of the
  resulting matrices contain contrasts which can be used for coding a
  factor with \code{n} levels.

  \code{contr.cum} gives a coding corresponding to successive differences
  between factor levels.

  \code{contr.diff} gives a coding that correspond to the cumulative sum
  of the value for each level. This is not meaningful in a model where the
  intercept is included, therefore \code{n} columns ia always returned.

  \code{contr.2nd} gives contrasts corresponding to 2nd order differences
  between factor levels. Returns a matrix with \code{n-2} columns.

  \code{contr.orth} gives a matrix with \code{n-2} columns, which are
  mutually orthogonal and orthogonal to the matrix \code{cbind(1,1:n)}
  }
\value{
  A matrix with \code{n} rows and \code{k} columns, with
  \code{k}=\code{n} for \code{contr.diff}
  \code{k}=\code{n-1} for \code{contr.cum}
  \code{k}=\code{n-2} for \code{contr.2nd} and  \code{contr.orth}.
  }
\author{Bendix Carstensen}
\seealso{
  \code{\link{contr.treatment}}
  }
\examples{
contr.cum(6)
contr.2nd(6)
contr.diff(6)
contr.orth(6)
  }
\keyword{design}
\keyword{models}
