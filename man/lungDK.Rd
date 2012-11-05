\name{lungDK}
\alias{lungDK}
\docType{data}
\title{Male lung cancer incidence in Denmark}
\description{
  Male lung cancer cases and population riks time in Denmark, for the
  period 1943--1992 in ages 40--89.
}
\usage{data(lungDK)}
\format{
  A data frame with 220 observations on the following 9 variables.
  \tabular{rl}{
    \code{A5}: \tab Left end point of the age interval, a numeric vector. \cr
    \code{P5}: \tab Left enpoint of the period interval, a numeric vector. \cr
    \code{C5}: \tab Left enpoint of the birth cohort interval, a numeric vector. \cr
    \code{up}: \tab Indicator of upper trianges of each age by period
      rectangle in the Lexis diagram. (\code{up=(P5-A5-C5)/5}). \cr
    \code{Ax}: \tab The mean age of diagnois (at risk) in the triangle. \cr
    \code{Px}: \tab The mean date of diagnosis (at risk) in the triangle. \cr
    \code{Cx}: \tab The mean date of birth in the triangle, a numeric vector. \cr
    \code{D}: \tab Number of diagnosed cases of male lung cancer. \cr
    \code{Y}: \tab Risk time in the male population, person-years. \cr
  }
}
\details{ Cases and person-years are tabulated by age and date of
  diagnosis (period) as well as date of birth (cohort) in 5-year
  classes. Each observation in the dataframe correponds to a triangle in
  a Lexis diagram. Triangles are classified by age and date of
  diagnosis, period of diagnosis and date of birth, all in 5-year
  groupings. 
}
\source{The Danish Cancer Registry and Statistics Denmark.
}
\references{
 For a more thorough exposition of statistical inference in the Lexis
 diagram, see: \url{http://staff.pubhealth.ku.dk/~bxc/APC/notes.pdf}
}
\examples{
data( lungDK )
# Draw a Lexis diagram and show the number of cases in it.
attach( lungDK )
Lexis.diagram( age=c(40,90), date=c(1943,1993), coh.grid=TRUE )
text( Px, Ax, paste( D ), cex=0.7 )
}
\keyword{datasets}
