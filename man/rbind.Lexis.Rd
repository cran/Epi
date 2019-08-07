\name{cbind.Lexis}
\alias{cbind.Lexis}
\alias{rbind.Lexis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combining a Lexis objects with data frames or other Lexis objects 
}
\description{A Lexis object may be combined side-by-side with
  data frames. Or several Lexis objects may stacked, possibly increasing
  the number of states and time scales.
}
\usage{
\method{cbind}{Lexis}(...)
\method{rbind}{Lexis}(...)
}
\arguments{
  \item{\dots}{For \code{cbind} a sequence of data frames or vectors of
    which exactly one has class \code{Lexis}. For \code{rbind} a sequence
    of Lexis objects, supposedly representing follow-up in the same
    population.}
  }
\details{
  Arguments to \code{rbind.Lexis} must all be \code{\link{Lexis}}
  objects; except for possible NULL objects. The timescales in the
  resulting object will be the union of all timescales present in all
  arguments. Values of timescales not present in a contributing Lexis
  object will be set to \code{NA}. The \code{breaks} for a given time
  scale will be \code{NULL} if the \code{breaks} of the same time scale
  from two contributing Lexis objects are different.

  The arguments to \code{cbind.Lexis} must consist of at most one Lexis
  object, so the method is intended for amending a Lexis object with
  extra columns without losing the Lexis-specific attributes.
}
\value{
  A \code{\link{Lexis}} object. \code{rbind} renders a \code{Lexis}
  object with timescales equal to the union of timescales in the
  arguments supplied. Values of a given timescale are set to \code{NA}
  for rows corresponding to supplied objects. \code{cbind} basically
  just adds columns to an existing Lexis object.
}
\author{
  Bendix Carstensen, \url{http://BendixCarstensen.com}
}
\seealso{
  \code{\link{subset.Lexis}}
}
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
xcoh <- cal.yr( xcoh, format="\%d/\%m/\%Y", wh=2:4 )
# See how it looks
xcoh
str( xcoh )

# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( per=entry ),
                exit = list( per=exit, age=exit-birth ),
         exit.status = fail,
                data = xcoh )
Lcoh
cbind( Lcoh, zz=3:5 )

# Lexis object wit time since entry time scale
Dcoh <- Lexis( entry = list( per=entry, tfe=0 ),
                exit = list( per=exit ),
         exit.status = fail,
                data = xcoh )
# A bit meningless to combie these two, really...
rbind( Dcoh, Lcoh )

# Split different places
sL <- splitLexis( Lcoh, time.scale="age", breaks=0:20*5 )
sD <- splitLexis( Dcoh, time.scale="tfe", breaks=0:50*2 )
sDL <- rbind( sD, sL )
}

\keyword{survival}
\keyword{manip}
