\name{splitLexis}
\alias{splitLexis}
\title{Split follow-up time in a Lexis object}
\description{
  The \code{splitLexis} function divides each row of a \code{Lexis}
  object into disjoint follow-up intervals according to the supplied
  break points.
}
\usage{
splitLexis(lex, breaks, time.scale, tol=.Machine$double.eps^0.5)
}
\arguments{
  \item{lex}{an object of class \code{Lexis}}
  \item{breaks}{a vector of break points}
  \item{time.scale}{the name or number of the time scale to be split}
  \item{tol}{numeric value >= 0. Intervals shorter than this value are
    dropped}
}
\value{
  An object of class \code{Lexis} with multiple rows for each row of
  the argument \code{lex}. Each row of the new \code{Lexis} object
  contains the part of the follow-up interval that falls inside one of
  the time bands defined by the break points.

  The variables representing the various time scales, are appropriately
  updated in the new \code{Lexis} object. The entry and exit status
  variables are also updated according to the rule that the entry status
  is retained until the end of follow-up. All other variables are
  considered to represent variables that are constant in time, and so
  are replicated across all rows having the same id value.
}
\note{
  The \code{splitLexis()} function divides follow-up time into intervals
  using breakpoints that are common to all rows of the \code{Lexis} object.
  To split a \code{Lexis} object by break points that are unique to each
  row, use the \code{cut.Lexis} function.
}
\author{Martyn Plummer}
\examples{
# A small bogus cohort
xcoh <- structure( list( id = c("A", "B", "C"),
                      birth = c("14/07/1952", "01/04/1954", "10/06/1987"),
                      entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                       exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                       fail = c(1, 0, 1) ),
                     .Names = c("id", "birth", "entry", "exit", "fail"),
                  row.names = c("1", "2", "3"),
                      class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh$bt <- cal.yr( xcoh$birth, format="\%d/\%m/\%Y" )
xcoh$en <- cal.yr( xcoh$entry, format="\%d/\%m/\%Y" )
xcoh$ex <- cal.yr( xcoh$exit , format="\%d/\%m/\%Y" )

# See how it looks
xcoh

# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( per=en ),
                exit = list( per=ex, age=ex-bt ),
         exit.status = fail,
                data = xcoh )

# Default plot of follow-up
plot( Lcoh )
# With a grid and deaths as endpoints
plot( Lcoh, grid=0:10*10, col="black" )
points( Lcoh, pch=c(NA,16)[Lcoh$lex.Xst+1] )
# With a lot of bells and whistles:
plot( Lcoh, grid=0:20*5, col="black", xaxs="i", yaxs="i",
      xlim=c(1960,2010), ylim=c(0,50), lwd=3, las=1 )
points( Lcoh, pch=c(NA,16)[Lcoh$lex.Xst+1], col="red", cex=1.5 )

# Split time along two time-axes
( x2 <- splitLexis( Lcoh, breaks = seq(1900,2000,5), time.scale="per") )
( x2 <- splitLexis( x2, breaks = seq(0,80,5), time.scale="age" ) ) 
str( x2 )

# Tabulate the cases and the person-years
tab( x2 )
tapply( status(x2,"exit")==1, list( timeBand(x2,"age","left"),
                                    timeBand(x2,"per","left") ), sum )
tapply( dur(x2),  list( timeBand(x2,"age","left"),
                        timeBand(x2,"per","left") ), sum )
}
\seealso{\code{\link{timeBand}}, \code{\link{cutLexis}}, \code{\link{tab.Lexis}}}
\keyword{manip}
