\name{cal.yr}
\alias{cal.yr}
\alias{as.Date.cal.yr}
\title{
  Functions to convert character, factor and various date objects into a number,
  and vice versa.
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
  cal.yr( x, format="\%Y-\%m-\%d", wh=NULL )
  as.Date.cal.yr( x, ... )
}
\arguments{
  \item{x}{A factor or character vector, representing a date in format
    \code{format}, or an object of class
    \code{\link{Date}},
    \code{\link{POSIXlt}},
    \code{\link{POSIXct}},
    \code{\link{date}},
    \code{dates} or
    \code{chron} (the latter two requires the \code{chron} package).
    If \code{x} is a data frame, all variables in the data-frame
    which are of one the classes mentioned are converted to class \code{cal.yr}.
    See arguemt \code{wh}, though.}
  \item{format}{Format of the date values if \code{x} is factor or character.
                If this argument is supplied and \code{x} is a datafame, all
                character variables are converted to class \code{cal.yr}.
                Factors in the dataframe will be ignored.}
  \item{wh}{Indices of the variables to convert if \code{x} is a data frame.
            Can be either a numerical or character vector.}
  \item{...}{Arguments passed on from other methods.}
  }
\value{
  \code{cal.yr} returns a numerical vector of the same length as
  \code{x}, of class \code{c("cal.yr","numeric")}. If \code{x} is a data frame
  a dataframe with some of the columns converted to class \code{"cal.yr"} is
  returned.

  \code{as.Date.cal.yr} returns a \code{\link{Date}} object.
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
 # Character vector of dates:
 birth <- c("14/07/1852","01/04/1954","10/06/1987","16/05/1990",
            "01/01/1996","01/01/1997","01/01/1998","01/01/1999")
 # Proper conversion to class "Date":
 birth.dat <- as.Date( birth, format="\%d/\%m/\%Y" )
 # Converson of character to class "cal.yr"
 bt.yr <- cal.yr( birth, format="\%d/\%m/\%Y" )
 # Back to class "Date":
 bt.dat <- as.Date( bt.yr )
 # Numerical calculation of days since 1.1.1970:
 days <- Days <- (bt.yr-1970)*365.25
 # Blunt assignment of class:
 class( Days ) <- "Date"
 # Then data.frame() to get readable output of results:
 data.frame( birth, birth.dat, bt.yr, bt.dat, days, Days, round(Days) )
}
\keyword{manip}
\keyword{chron}
