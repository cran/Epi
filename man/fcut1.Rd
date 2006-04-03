\name{fcut1}
\alias{fcut1}
\title{
  Cut follow-up time at a failure time.
}
\description{
  This function cuts the follow-up time at a failure allowing a person
  to stay at risk after the failure. It is aimed at processing of
  recurrent events. Failure times outside the interval
  (\code{enter},\code{exit}) are ignored.
}
\usage{
fcut1( enter, exit, fail, dof, fail.value = 1, 
         data = data.frame(enter, exit, fail, dof),
       Expand = 1:nrow(data) )
}
\arguments{
  \item{enter}{Date of entry into the study. Numerical vector.}
  \item{exit}{Date of exit from the study. Numerical vector.}
  \item{fail}{Failure indicator for the exit date.}
  \item{dof}{Date of failure. Numerical vector. Missing for persons who
    do not have an event.}
  \item{fail.value}{Value for failure indicator \code{Fail} at the date
    \code{dof}.}
  \item{data}{Dataframe of variables to be carried unchanged to the
    output.}
  \item{Expand}{Variable identifying original records.}
}
\value{
  A dataframe with the same columns as in \code{data}, preceded by the columns:
  \item{Expand}{Identification of the rows from the input dataframe.}
  \item{Enter}{Entry date for the interval.}
  \item{Exit}{Exit date for the interval.}
  \item{Fail}{Failure indicator for the interval.}
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \email{bxc@steno.dk}, \url{www.biostat.ku.dk/~bxc}
}
\seealso{
  \code{\link{Lexis}},
  \code{\link{isec}},
  \code{\link{fcut1}},
  \code{\link{fcut}},
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
# Some failure indicators and failure times
fail <- sample( 0:1, 15, replace=TRUE, prob=c(0.7,0.3) )
dof <- sample( c(one,two), 15 )
# So what have we got
data.frame( doe, dox, fail, dof )
# Cut follow-up at dof
fcut1( doe, dox, fail, dof )
}
\keyword{manip}
\keyword{datagen}
