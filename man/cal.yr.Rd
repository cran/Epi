\name{Calendar year}
\alias{cal.yr}
\alias{dmy.cal.yr}
\alias{mdy.cal.yr}
\alias{ymd.cal.yr}
\alias{weekdays.cal.yr}
\alias{weeks}
\alias{weeks.cal.yr}
\alias{months.cal.yr}
\alias{print.cal.yr}
\title{
  Functions representing dates as fractional calendar years
  }
\description{
  Dates are converted to a numerical value, so that 1 January
  1970 is converted to 1970.0. The result is given a class
  attribute, \code{cal.yr}. The function is intended for epidemiological
  use, and the dates are converted to a calendar time scale that has 4
  years of equal length between 1.1.1970 and 1.1.1974 etc. The printing of
  \code{cal.yr} objects in human readable form may therefore be as much
  as 3 days wrong.
  }
\usage{
  cal.yr(x, format)
  dmy.cal.yr( d, m, y )
  mdy.cal.yr( m, d, y )
  ymd.cal.yr( y, m, d )
  weekdays.cal.yr( x, abbreviate=FALSE )
  weeks( x )
  weeks.cal.yr( x )
  months.cal.yr( x, abbreviate=FALSE )
  print.cal.yr( x, format="\%d/\%m/\%Y", ... )
}
\arguments{
  \item{x}{A character vector, representing a date in format
    \code{format}, or an object of class \code{\link{POSIXt}} or
    \code{\link{Date}}.} 
  \item{d, m, y}{Numerical vectors representing, date month and year.}
  \item{format}{Optional character values giving the format for the date
    if given as character. See the help page for \code{\link{strptime}} for
    admissible values.}
  \item{abbreviate}{Should names be abbreviated?}
  \item{...}{Ignored, but present for compatibility with print methods}
  }
\value{
  \code{cal.yr} returns a numerical vector of the same length as
  \code{x}. The associated print method, \code{print.cal.yr} prints the
  dates in human readable format. To see the dates as fractions of
  years, use \code{as.numeric(cy)}.

  \code{dmy.cal.yr}, \code{mdy.cal.yr} and \code{ymd.cal.yr} are
  convenience functions for converting separate numerical variables
  with date \code{d}, month \code{m} and year \code{y} into dates of
  class \code{cal.yr}. These functions also accept a single character
  variable as input, with two digits for day and month and four for year
  and any delimiters.

  \code{weekdays.cal.yr}, \code{weeks.cal.yr} and \code{months.cal.yr},
  extract weekday, week number and month from calendar year
  objects. \code{weekdays} and \code{months} return factors with 7 and
  12 levels respectively, regardless of the actually occurring number of
  levels. \code{weeks} return numeric values in the range 1--52.
}
\author{
  Bendix Carstensen, Steno Diabetes Center \& Dept. of Biostatistics,
  University of Copenhagen, \url{http://www.biostat.ku.dk/~bxc}
} 
\seealso{
  \code{\link{DateTimeClasses}},
  \code{\link{Date}}
}
\examples{
 birth <- c("14/07/1952", "01/04/1954", "10/06/1987")
 ( bt.yr <- cal.yr( birth, format="\%d/\%m/\%Y" ) )
 as.numeric( bt.yr )
 weeks( bt.yr )
 months( bt.yr )
 ( bt.y <- dmy.cal.yr( birth ) )
}
\keyword{manip}
\keyword{chron}
