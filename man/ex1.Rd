\name{ex1}
\alias{ex1}
\title{
  Split follow-up time along a timescale
}
\description{
  Splits follow-up time at prespecified points of follow-up.
}
\usage{
ex1( enter, exit, fail, origin = 0, scale = 1, breaks,
     data = data.frame(enter, exit, fail),
     Expand = 1:nrow(data) )
}
\arguments{
  \item{enter}{Date of entry einto the study (start of follow-up). Numeric.}
  \item{exit}{Date of termination of follow-up. Numeric.}
  \item{fail}{Status at exit from the study.}
  \item{origin}{Origin of the timescale to split on. Specified on the
    input timescale, i.e. that of \code{enter} and \code{exit}.}
  \item{scale}{Scaling between input and analysis timescale.}
  \item{breaks}{Breakpoints on the analysis timescale. Follow-up before
    \code{min(breaks)} and and ater \code{max(breks)} is discarded.}
  \item{data}{Dataframe of variables to carry over to the output.}
  \item{Expand}{Variable identifying original records.}
}
\details{
  If \code{entry} and \code{exit} are given in days (for example as
  \code{Date} variables, and we want follow-up cut at 5-year age
  intervals, the we should choose \code{origin} equal to bithdate, scale
  equal to 365.25 and breaks as \code{seq(0,100,5)}. Thus the input
  timescale is calendar tiem measured in days, and output timescale is
  age measured in years. 
}
\value{
  A dataframe with one row per follow-up interval, and variables as in
\code{data}, preceded by the variables:
  \item{Expand}{Identification of the rows from the input dataframe.}
  \item{Enter}{Entry date for the interval.}
  \item{Exit}{Exit date for the interval.}
  \item{Fail}{Failure indicator for end of the current interval.}  
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
# Split follow-up:
ex1( doe, dox, fail, breaks=0:10 )
}
\keyword{manip}
\keyword{datagen}
