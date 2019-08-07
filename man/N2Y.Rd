\name{N2Y}
\alias{N2Y}
\title{
Create risk time ("Person-Years") in Lexis triangles from population count data.
}
\description{
Data on population size at equidistant dates and age-classes are
used to estimate person-time at risk in Lexis-triangles, i.e. classes
classified by age, period AND cohort (date of birth). Only works for
data where age-classes have the same width as the period-intervals.
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
\value{A data frame with variables \code{A}, \code{P} and \code{Y},
  representing the mean age and period in the Lexis triangles and the
  person-time in them, respectively. The person-time is in units of the
  distance between population count dates.

  If \code{res.dfr=FALSE} a three-way table classified by the left end
  point of the age-classes and the periods and a factor \code{wh} taking
  the values \code{up} and \code{lo} corresponding to upper (early
  cohort) and lower (late cohort) Lexis triangles.
}
\details{The calculation of the risk time from the population figures is
  done as described in: B. Carstensen: Age-Period-Cohort models for the
  Lexis diagram. Statistics in Medicine, 26: 3018-3045, 2007.

  The number of periods in the result is one less than the number
  of dates (\code{nP=length(table(P))}) in the input, so the number of
  distinct values is \code{2*(nP-1)}, because the \code{P} in the output
  is coded differently for upper and lower Lexis triangles.

  The number of age-classes is the same as in the input.  In the paper
  "Age-Period-Cohort models for the Lexis diagram" I suggest that the
  risk time in the lower triangles in the first age-class and in the
  upper triangles in the last age-class are computed so that the total
  risk time in the age-class corresponds to the average of the two
  population figures for the age-class at either end of the period. This
  is the method used.
}
\references{
  B. Carstensen: Age-Period-Cohort models for the
  Lexis diagram. Statistics in Medicine, 26: 3018-3045, 2007.
}
\author{
  Bendix Carstensen, \url{http://BendixCarstensen.com}
       }
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

# Calculation of PY for persons born 1970 in 1972
( N.1.1972 <- subset( Nx, A==1 & P==1972)$N )
( N.2.1973 <- subset( Nx, A==2 & P==1973)$N )
N.1.1972/3 + N.2.1973/6
N.1.1972/6 + N.2.1973/3
# These number can be found in the following plot:

# Blue numbers are population size at 1 January
# Red numbers are the computed person-years in Lexis triangles:
Lexis.diagram(age=c(0,5), date=c(1970,1975), int=1, coh.grid=TRUE )
with( Nx, text(P,A+0.5,paste(N),srt=90,col="blue") )
with( Nt, text(P,A,formatC(Y,format="f",digits=1),col="red") )
text( 1970.5, 2, "Population count 1 January", srt=90, col="blue")
text( 1974.5, 2, "Person-\nyears", col="red")
}
\keyword{Data}
