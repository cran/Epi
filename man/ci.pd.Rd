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
ci.pd(aa, bb=NULL, cc=NULL, dd=NULL,
     method = "Nc",
      alpha = 0.05, conf.level=0.95,
     digits = 3,
      print = TRUE,
detail.labs = FALSE )
}
\arguments{
  \item{aa}{Numeric vector of successes in sample 1. Can also be a
    matrix or array (see details).} 
  \item{bb}{Successes in sample 2.}
  \item{cc}{Failures in sample 1.}
  \item{dd}{Failures in sample 2.}
  \item{method}{Method to use for calculation of confidence interval, see
    "Details".}
  \item{alpha}{Significance level}
  \item{conf.level}{Confidence level}
  \item{print}{Should an account of the two by two table be printed.}
  \item{digits}{How many digits should the result be rounded to if printed.}
  \item{detail.labs}{Should the computing of probability differences be
                     reported in the labels.}
}
\details{
  Implements method 10 from Newcombe(1998) (method="Nc") or from
  Agresti & Caffo(2000) (method="AC").

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
  
  A Agresti & B Caffo: Simple and effective confidence intervals for
  proportions and differences of proportions result from adding two
  successes and two failures. The American Statistician,
  54(4), pp. 280-288, 2000.
}
\author{
  Bendix Carstensen, Esa L‰‰r‰.
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
ci.pd( A )
ci.pd( A, detail.labs=TRUE, digits=3 )
}
\keyword{distribution}
\keyword{htest}
