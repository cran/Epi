\name{in.span}
\alias{in.span}
\alias{inSpan}
\alias{id.span}
\alias{idSpan}
\alias{thinCol}
\title{
Is \code{x} in the column span of matrix \code{A} and what columns are
linearly dependent?
}
\description{
The function \code{in.span} checks if the vector \code{x} (or
columns of the matrix \code{x}) is in the column span of the matrix
\code{A}. If desired, it returns the coefficient matrix \code{B} so that
\code{AB=x}. The function \code{thinCol} removes linearly dependent
columns an returns a matrix of full rank.
}
\usage{
in.span( A,
         x,
      coef = FALSE,
       tol = 1e-08 )
inSpan( A, x, coef=FALSE, tol=1e-08 )
id.span( A, B, tol=1e-08 )
idSpan( A, B, tol=1e-08 )
thinCol( A, tol = 1e-06, col.num = FALSE )
}
\details{\code{\link{thinCol}} is mainly a workhorse in
  \code{\link{detrend}}, but made available because of its general
  usefulness. 

  \code{in.span} and \code{inSpan} are just different names for the same
  to accommodate different naming schools.

  \code{in.span} (\code{inSpan}) is handy in checking whether different
  parametrizations of a model are identical in the sense of spanning the
  same linear space. Equivalent to checking whether fitted values under
  different parametrizations are identical, but has the further use of
  checking if subspaces of models are equivalent. The function
  simply checks if the regression of (columns of) \code{x} on the
  columns of \code{A} produces residuals that are all 0.

  \code{id.span} (equivalent to \code{idSpan}) checks whether two
  matrices have the same column span.  
}
\arguments{
\item{A}{A matrix.}
\item{B}{A matrix.}
\item{x}{A vector or matrix. \code{length(x)} (or \code{nrow(x)}) must
  be equal to \code{nrow(A)}.} 
\item{coef}{Logical. Should the coefficient matrix (\code{k}) be
  returned, so that \code{Ak=x}?}
\item{tol}{Tolerance for identity of matrices in check
  (\code{in.span}) or QR decomposition (\code{thinCol})}
\item{col.num}{Logical. Should the positions of dependent columns be
  returned instead of the full-rank matrix?}
}
\value{\code{in.span} returns a logical: is \code{x} is in
  \code{span(A)}? If \code{coef=TRUE} it returns a matrix \code{k} so
  that \code{Ak=x}. \code{k} is not necessarily unique (A may not have
  full rank).

  \code{id.span} returns a logical: is \code{span(A)} the same as
  \code{span(B)}? 

  \code{thinCol} returns a matrix of full rank, formed from \code{A} by
  deleting columns linearly dependent on other. If \code{col.num=TRUE}
  (one possible set of) positions of columns forming a full rank basis
  for the column space of \code{A} is returned.
}
\author{Bendix Carstensen, \url{BendixCarstensen.com} with essential
  help from Lars Jorge Diaz and Peter Dalgaard.
}
\seealso{\code{\link{det}}
}
\examples{
# Matrices and vectors, x in span(A), z (hopefully) not
A <- matrix(round(rnorm(15)*20),5,3)
B <- matrix(round(rnorm(15)*20),5,3)
B <- cbind( B, B\%*\%c(3,4,2) )
x <- A \%*\% c(3,4,2)
z <- 5:9

# how they look
data.frame( A=A, x=x, z=z, B=B )

# vectors in span(A)?
in.span(A,x)
in.span(x,A)
in.span(A,x,coef=TRUE)

in.span(A,z)
in.span(A,z,coef=TRUE)

# Do matrices span the same space ?
in.span( A, B )
in.span( B, A )

# B is not in span of a subspace of B columns, but vice versa
( M <- matrix( rnorm(8)*7, 4, 2 ) )
in.span( B\%*\%M, B )
in.span( B, B\%*\%M )
id.span( B, B\%*\%M )

# But not unique for singular matrices:
( xx <- in.span( B, B\%*\%M, coef=TRUE ) )
cbind( B\%*\%M, B\%*\%xx )
cbind( xx, M )

# Easier for full rank matrices:
( K <- matrix( rnorm(9)*7, 3, 3 ) )
in.span( A\%*\%K, A )
in.span( A, A\%*\%K )
id.span( A, A\%*\%K )
in.span( A, A\%*\%K, coef=TRUE )
}
\keyword{math}
