\name{isec}
\alias{isec}
\title{
  Determine the intersection between follow-up intervals and a fixed
  interval.
}
\description{
For a given piece of follow-up for ach person, this function determines
the part of the follow-up which is inside a pre-specified interval. The
function may occasionally be useful but is mainly included as a tool to
be used in the function \code{\link{Lexis}}.
}
\usage{
  isec( enter, exit, fail = 0, int, cens.value = 0,
        Expand = 1:length(enter))
}
\arguments{
  \item{enter}{Numereical vector of entry times.}
  \item{exit}{Numerical vector of exit times.}
  \item{fail}{Vector of exit status, i.e. status at time \code{exit}}
  \item{int}{The fixed interval. Numerical vector of length 2.}
  \item{cens.value}{Censoring value. Numerical. If the person survives the fixed
    interval, this will be the value of the status indicator, \code{Fail}}
  \item{Expand}{Person id to be carried to the output}
}
\value{
  A matrix with columns \code{Expand} - person identification,
  \code{Enter} - date of entry into the interval, \code{Exit} - date of exit
  from the interval and \code{Fail} - status indicator. 
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \email{bxc@steno.dk}, \url{www.biostat.ku.dk/~bxc}
}
\seealso{
  \code{\link{Lexis}},
  \code{\link{icut}},
  \code{\link{fcut1}},
  \code{\link{fcut}},
  \code{\link{ex1}}
}
\examples{
one <- round( runif( 15, 0, 15 ), 1 )
two <- round( runif( 15, 0, 15 ), 1 )
doe <- pmin( one, two )
dox <- pmax( one, two )
# Goofy data rows to test possibly odd behaviour
doe[1:3] <- dox[1:3] <- 8
dox[2] <- 6
dox[3] <- 7.5
# Some failure indicators
fail <- sample( 0:1, 15, replace=TRUE, prob=c(0.7,0.3) )
# So what have we got?
data.frame( doe, dox, fail )
# Find intersection with interval (4,8)
isec( doe, dox, fail, int=c(4,8) )
# See how it compares to original data
merge( data.frame( Expand=1:15, doe, dox, fail ),
       data.frame( isec( doe, dox, fail, int=c(4,8) ) ), all=TRUE )
  }
}
\keyword{manip}
\keyword{datagen}
