\name{icut}
\alias{icut}
\title{
  Function to cut the follow-up in cohort at a point in time.
}
\description{
  The follow up time from \code{enter} to \code{exit} is classified as
  to wheter it is before \code{cut} (\code{Time}=0) or after
  \code{Time}=1). If \code{cut} is between \code{enter} and \code{exit},
  the follow-up is split in two intervals, the first gets the value
  \code{cens} for the status.
}
\usage{
icut( enter, exit, cut, fail = 0, cens.value = 0,
      data = data.frame(enter, exit, fail, cut),
      Expand = 1:nrow( data ),
      na.cut = Inf )
}
\arguments{
  \item{enter}{Date of entry. Numerical vector.}
  \item{exit}{Date of exit. Numerical vector.}
  \item{fail}{Indicator if exit status.}
  \item{cens.value}{Value to use for censoring status.}
  \item{cut}{Date where to cut follow-up. Numerical vector.}
  \item{data}{Dataframe of variables to carry over to the output dataframe.}
  \item{Expand}{Variable identifying original records.}
  \item{na.cut}{What value should be assigned to missing values of the
    cutpoint. Defaults to \code{Inf}, so the inetrmediate event is
    considered not to have occcurred. If set to \code{-Inf}, all persons
    with missing \code{cut} are considered to have had an intermediate
    event. If set to \code{NA} records with missing \code{cut} are
    omitted from the result.} 
}
\details{
  The purpose of this function is to divide follow-up into pre- and post
  some intermediate event like recurrence of disease, thus enabling
  Follow-up for persons with a recurrence date (\code{cut}) will be
  split in two, with indication (in \code{Time}) of what is pre and what
  is post recurrence. This is typically what precedes a survival
  analysis where recurrence is modelled as a time-dependent variable.
}
\value{
  A data frame with one row per interval of follow up and columns given
  in the \code{data} argument, preceded by the columns:
  \item{Expand}{Identification of the rows from the input dataframe.}
  \item{Enter}{Entry date for the interval.}
  \item{Exit}{Exit date for the interval.}
  \item{Fail}{Failure indicator for the interval.}
  \item{Time}{Indicator variable for intervals after \code{cut}.}
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
# So what have we got
data.frame( doe, dox, fail )
# Cut follow-up at 5
icut( doe, dox, fail, cut=5 )
}
\keyword{manip}
\keyword{datagen}
