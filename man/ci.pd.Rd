\name{ci.pd}
\alias{ci.pd}
\title{
  Compute confidence limits for a difference of two independent proportions.
}
\description{
  The usual formula for the c.i. of at difference of proportions is
  inaccurate. Newcombe has compared 11 methods and method 10 in his
  paper looks like a winner. It is implemented here.
}
\usage{
ci.pd(aa, bb, cc, dd, alpha = 0.05, print = TRUE)
}
\arguments{
  \item{aa}{Numeric vector of successes in sample 1. Can also be a
    matrix (see details).} 
  \item{bb}{Successes in sample 2.}
  \item{cc}{Failures in sample 1.}
  \item{dd}{Failures in sample 2.}
  \item{alpha}{Significance level}
  \item{print}{Should an account of the two by two table be
    printed. Ignored if more than difference is computed, i.e. if
    \code{aa}, \code{bb}, \code{cc} and \code{dd} are vectors or if \code{aa}
    is a 3-dimensional matrix.} 
}
\details{
  \code{aa}, \code{bb}, \code{cc} and \code{dd} can be vectors.
  If \code{aa} is a matrix, the elements \code{[1:2,1:2]} are used, with
  successes \code{aa[,1:2]}. If \code{aa} is a three-way table or array,
  the elements \code{aa[1:2,1:2,]} are used.
}
\value{
  A matrix with three columns: probability difference, lower and upper
  limit. The number of rows equals the length of the vectors  \code{aa},
  \code{bb}, \code{cc} and \code{dd} or, if \code{aa} is a 3-way matrix,
  \code{dim(aa)[3]}.   
}
\references{
  RG Newcombe: Interval estimation for the difference between
  independent proportions. Comparison of eleven methods. Statistics in
  Medicine, 17, pp. 873-890, 1998.
}
\author{
  Bendix Carstensen,
  \url{http://www.biostat.ku.dk/~bxc}
}
\seealso{
  \code{\link{twoby2}}, \code{\link{binom.test}}
}
\examples{
( a <- matrix( sample( 10:40, 4 ), 2, 2 ) )
ci.pd( a )
twoby2( t(a) )
prop.test( t(a) )
( A <- array( sample( 10:40, 20 ), dim=c(2,2,5) ) )
print( ci.pd( A ) )
}
\keyword{distribution}
\keyword{htest}
