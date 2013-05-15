\name{N2Y}
\alias{N2Y}
\title{
Create risk time (Person-Years) in Lexis triangles from population data.
}
\description{
Data on population size at equidistant dates and age-classes are
used to estimate person-time at risk in Lexis-triangles, i.e. classes
classified by age, period AND cohort.
}
\usage{
N2Y( A, P, N,
     data = NULL,
     return.dfr = TRUE)
}
\arguments{
  \item{A}{Name of the age-variable, which should be numeric,
    corresponding to the left endpoints of the age intervals.}
  \item{P}{Name of the period-variable, which should be numeric,
    corresponding to the date of population count.}
  \item{N}{The population size at date \code{P} in age class \code{A}.}
  \item{data}{A data frame in which arguments are interpreted.}
  \item{return.dfr}{Logical. Should the results be returned as a data frame
    (default \code{TRUE}) or as a table.}
}
\details{The calculation of the risk time from the population figures is
  done as described in: B. Carstensen: Age-Period-Cohort models for the
  Lexis diagram. Statistics in Medicine, 26: 3018-3045, 2007.
}
\value{A data frame with variables \code{A}, \code{P} and \code{Y},
  representing the mean age and period in the Lexis triangles and the
  person-time in them.

  If \code{res.dfr=FALSE} a three-way table classified by the left end
  point of the age-classes and the periods and a factor \code{wh} taking
  the values \code{up} and \code{lo} corresponding to upper (early
  cohort) and lower (late cohort) Lexis triangles.
}
\references{
  B. Carstensen: Age-Period-Cohort models for the
  Lexis diagram. Statistics in Medicine, 26: 3018-3045, 2007.
}
\author{  Bendix Carstensen,
  \url{BendixCarstensen.com} }
\seealso{
\code{\link{splitLexis}}, \code{\link{apc.fit}}
}
\examples{
# Danish population at 1 Jan each year by sex and age
data( N.dk )
# An illustrative subset
( Nx <- subset( N.dk, sex==1 & A<5 & P<1975 ) )
# Show the data in tabular form
xtabs( N ~ A + P, data=Nx )
# Lexis triangles as data frame
Nt <- N2Y( data=Nx, return.dfr=TRUE )
xtabs( Y ~ round(A,2) + round(P,2), data=Nt )
# Lexis triangles as a 3-dim array
ftable( N2Y( data=Nx, return.dfr=FALSE ) )
}
\keyword{Data}
