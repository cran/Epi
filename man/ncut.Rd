\name{ncut}
\alias{ncut}
\title{
Function to group a variable in intervals.}
\description{
  Cuts a continuous variable in intervals. As opposed to \code{cut}
  which returns a factor, \code{ncut} returns a numeric variable.
}
\usage{
ncut(x, breaks, type="left" )
}
\arguments{
  \item{x}{A numerical vector.}
  \item{breaks}{Vector of breakpoints. \code{NA} will results for values
    below \code{min(x)} if \code{type="left"}, for values
    above \code{max(x)} if \code{type="right"} and for values
    outside \code{range(x)} if \code{type="mid"}}
  \item{type}{Character: one of \code{c("left","right","mid")},
    indicating whether the left, right or midpoint of the intervals
    defined in breaks is returned.}
}
\details{
The function uses the base function \code{findInterval}.
}
\value{
A numerical vector of the same length as \code{x}.
}
\author{
  Bendix Carstensen, Steno Diabetes Center, \email{bxc@steno.dk},
  \url{http://www.biostat.ku.dk/~bxc/}, with essential input
  from Martyn Plummer, IARC.
}
\seealso{
  \code{\link{cut}}, \code{\link{findInterval}}
}
\examples{
br <- c(-2,0,1,2.5)
x <- c( rnorm( 10 ), br, -3, 3 )
cbind( x, l=ncut( x, breaks=br, type="l" ),
          m=ncut( x, breaks=br, type="m" ),
          r=ncut( x, breaks=br, type="r" ) )[order(x),]
x <- rnorm( 200 )
plot( x, ncut( x, breaks=br, type="l" ), pch=16, col="blue", ylim=range(x) )
abline( 0, 1 )
abline( v=br )
points( x, ncut( x, breaks=br, type="r" ), pch=16, col="red" )
points( x, ncut( x, breaks=br, type="m" ), pch=16, col="green" )
}
\keyword{manip}

