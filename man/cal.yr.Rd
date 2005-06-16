\name{cal.yr}
\alias{cal.yr}
\title{
  Function to convert character, factor and various date objects into a number.
  }
\description{
  Dates are converted to a numerical value, giving the calendar year as
  a fractional number. 1 January 1970 is converted to 1970.0, and other
  dates are converted by assuming that years are all 365.25 days long,
  so inaccuracies may arise, for example, 1 Jan 2000 is converted to
  1999.999. Differences between converted values will be 1/365.25 of the 
  difference between corresponding \code{\link{Date}} objects.
  }
\usage{
  cal.yr(x, format)
}
\arguments{
  \item{x}{A factor or character vector, representing a date in format
    \code{format}, or an object of class
    \code{\link{Date}},
    \code{\link{POSIXlt}},
    \code{\link{POSIXct}},
    \code{\link{date}},
    \code{dates} or
    \code{chron} (the latter two requires the \code{chron} package).}
  \item{format}{Format of the date values if \code{x} is factor or character}
  }
\value{
  \code{cal.yr} returns a numerical vector of the same length as
  \code{x}.
}
\author{
  Bendix Carstensen, Steno Diabetes Center \& Dept. of Biostatistics,
  University of Copenhagen, \email{bxc@steno.dk},
  \url{http://www.pubhealth.ku.dk/~bxc}
} 
\seealso{
  \code{\link{DateTimeClasses}},
  \code{\link{Date}}
}
\examples{
 birth <- c("14/07/1952","01/04/1954","10/06/1987","16/05/1990",
            "01/01/1996","01/01/1997","01/01/1998","01/01/1999")
 ( bt.yr <- cal.yr( birth, format="\%d/\%m/\%Y" ) )
 # To get readable output:
 data.frame( birth, bt.yr )
}
\keyword{manip}
\keyword{chron}
