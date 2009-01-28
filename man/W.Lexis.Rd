\name{W.Lexis}
\alias{W.Lexis}
\title{Split follow-up time in cohort studies.}
\description{
  This is an old version of time-splitting. Consult the new version,
  \code{\link{Lexis}}.
  
  For cohort input data the follow-up time is chopped in pieces along
  several time scales, and a dataframe of follow-up intervals is
  returned. Entry and exit times are assumed to be in the same timescale
  (the input time scale).
}
\usage{
W.Lexis( entry = 0,
          exit,
          fail,
        origin = 0,
         scale = 1,
        breaks,
       include = NULL,
          data = NULL )
}
\arguments{
  \item{entry}{Date of entry on the input timescale. Numerical variable.}
  \item{exit}{Date of exit on the input timescale. Numerical variable.}
  \item{fail}{Failure indicator.}
  \item{origin}{Origin of the output timescale(s) on the input
    timescale. If for example the input timescale is calendar time and
    the output timescale is (current) age, the the origin is date of
    birth. If more than one timescale is used for splitting time
    this is a list. Elements of the list must be named and must have the
    same names as those in \code{scale} and \code{breaks}.} 
  \item{scale}{Scale of the output timescale(s) relative to the input
    timescale. Elements of the list must be named and have the same
    names as those in \code{origin} and \code{breaks}.} 
  \item{breaks}{Points on the output scale where the follow-up is
    cut. If more than one timescale is used for splitting time this is a
    list. Elements of the list must be named and must have the same
    names as those in \code{origin} and \code{scale}.} 
  \item{include}{List of variables to carry unchanged from the original
    dataframe to the output dataframe.} 
  \item{data}{Dataframe in which to interpret the arguments.}
}
\details{
  The \code{data} is assumed to be a dataframe describing the follow-up
  of a cohort, giving entry and exit time (on the input timescale) for each
  individual as well as the exit status (failure status,
  \code{fail}). The purpose of the function is to split each
  individual's follow-up time along a number of timescales for example
  age, calendar time, time since entry etc.
  Any follow-up time before the first break
  on any timescale or after the last break on any of these timescales
  (the output timescales) is discarded.

  NOTE: If a person has his/her exit before the first break or his entry
  after the last break on any of the timescales the function will crash.
}
\value{
  A dataframe with one row per follow-up interval. The data frame includes
  all variables given in the argument \code{include}, replicated across
  all intervals from the same individual.  In addition, the dataframe has
  the following variables:
  \item{Expand}{A numerical vector with values in
    \code{1:nrows(data)}, pointing at the rows of the input data frame
    that is expanded to the output intervals.}
  \item{Entry}{Date of entry for each interval. On the input time scale.}
  \item{Exit}{Date of exit for each interval. On the input time scale.}
  \item{Fail}{Exit status for each interval. Coded 0 for censoring, for
    the last follow-up interval for each person it takes the value of
    \code{fail}.} 
  \item{Time}{If \code{origin}, \code{scale} or \code{breaks} were given
    as vectors this gives the left endpoints of the intervals
    on the output scale.

    If \code{origin}, \code{scale} or
    \code{breaks}, were given as lists, there is no variable \code{Time}
    in the dataframe, instead variables with the same names as
    the list elements of these will be in the dataframe. The variables
    have values corresponding to the left endpoints of the intervals on
    the respective output time scales.}
}
\references{
  This function has approximately the same functionality as
  \code{stsplit} in Stata and the SAS-macro \code{\%Lexis}
  (\url{http://www.biostat.ku.dk/~bxc/Lexis/Lexis.sas}). It has been
  attempted to keep argument names similar between the three functions.
}  
\author{
  David Clayton, approx. 2000.
  Small modifications by Bendix Carstensen.
}
\seealso{
  \code{\link{Lexis}},\code{\link{Lexis.diagram}}
}
\examples{
# A small bogus cohort
#
xcoh <- structure( list( id = c("A", "B", "C"),
                      birth = c("14/07/1952", "01/04/1954", "10/06/1987"),
                      entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                       exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                       fail = c(1, 0, 1) ),
                     .Names = c("id", "birth", "entry", "exit", "fail"),
                  row.names = c("1", "2", "3"),
                      class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
#
xcoh$bt <- cal.yr( xcoh$birth, format="\%d/\%m/\%Y" )
xcoh$en <- cal.yr( xcoh$entry, format="\%d/\%m/\%Y" )
xcoh$ex <- cal.yr( xcoh$exit , format="\%d/\%m/\%Y" )

# See how it looks
#
xcoh 

# Split time along one time-axis
#
W.Lexis( entry = en,
          exit = ex,
          fail = fail,
         scale = 1,
        origin = bt,
        breaks = seq( 5, 40, 5 ),
       include = list( bt, en, ex, id ),
          data = xcoh )

# Split time along two time-axes
#
( x2 <- 
W.Lexis( entry = en,
          exit = ex,
          fail = fail,
         scale = 1,
        origin = list( per=0,                 age=bt          ),
        breaks = list( per=seq(1900,2000,10), age=seq(0,80,5) ),
       include = list( bt, en, ex, id ),
          data = xcoh ) )

# Tabulate the cases and the person-years
#
tapply( x2$Fail, list( x2$age, x2$per ), sum )
tapply( x2$Exit - x2$Entry, list( x2$age, x2$per ), sum )
}
\keyword{manip}
