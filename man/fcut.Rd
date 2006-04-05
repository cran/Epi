\name{fcut}
\alias{fcut}
\title{
  Cuts follow-up time at multiple failure times.
}
\description{
  This function cuts the follow-up time at multiple failure times
  allowing a person to stay at risk between and after the laset
  failure. It is aimed at processing of recurrent events. Failure times
  outside the interval (\code{enter},\code{exit}) are ignored.
}
\usage{
fcut( enter, exit, dof, fail = 0,
      data = data.frame(enter, exit),
      Expand = 1:nrow( data ))
}
\arguments{
  \item{enter}{Date of entry into the study. Numerical vector.}
  \item{exit}{Date of exit from the study. Numerical vector.}
  \item{fail}{Failure indicator for the exit date.}
  \item{dof}{Failure time(s). For multiple failures per individual,
    \code{dof} must be a list.}
  \item{data}{Dataframe of variables to be carried over to the output.}
  \item{Expand}{Variable identifying original records.}
}
\value{
  A dataframe with the same variables as in \code{data} preceded by the
  variables:
  \item{Expand}{Identification of the rows from the input dataframe.}
  \item{Enter}{Entry date for the interval.}
  \item{Exit}{Exit date for the interval.}
  \item{Fail}{Failure indicator for end of the current interval.}
  \item{n.Fail}{Number of failures prior to the start of the current
    interval. Counts all failures given in the list \code{dof},
    including those prior to \code{enter}.}
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \email{bxc@steno.dk}, \url{www.biostat.ku.dk/~bxc}
}
\seealso{
  \code{\link{Lexis}},
  \code{\link{isec}},
  \code{\link{icut}},
  \code{\link{fcut1}},
  \code{\link{ex1}}
}
\examples{
one <- round( runif( 15, 0, 10 ), 1 )
two <- round( runif( 15, 0, 10 ), 1 )
doe <- pmin( one, two )
dox <- pmax( one, two )
# Goofy data rows to test possibly odd behaviour
doe[1:3] <- dox[1:3] <- 8
dox[2] <- 6
dox[3] <- 7.5
# Some failure indicators
fail <- sample( 0:1, 15, replace=TRUE, prob=c(0.7,0.3) )
# Failure times in a list
dof <- sample( c(one,two), 15 )
l.dof <- list( f1=sample( c(one,two), 15 ),
               f2=sample( c(one,two), 15 ),
               f3=sample( c(one,two),15 ) )
# The same, but with events prior to entry removed
lx.dof <- lapply( l.dof, FUN=function(x){ x[x<doe] <- NA ; x } )
# So what have we got
data.frame( doe, dox, fail, l.dof, lx.dof )
# Cut follow-up at event times
fcut( doe, dox, lx.dof, fail, data=data.frame( doe, dox, lx.dof ) )
}
\keyword{manip}
\keyword{datagen}
