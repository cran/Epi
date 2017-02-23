\name{DMepi}
\alias{DMepi}
\docType{data}
\title{Epidmiological rates for diabetes in Denmark 1996--2015}
\description{Register based counts and person-uears for incidece of
  diabetes and mortality with and without diabetes.
}
\usage{data("DMepi")}
\format{
  A data frame with 4000 observations on the following 8 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{M} \code{F}}
    \item{\code{A}}{Age glass 0 -- 99}
    \item{\code{P}}{Calendar year, 1996-2015}
    \item{\code{X}}{Number of new diagnoses of diabetes}
    \item{\code{D.nD}}{Number of deaths among persons without diabetes}	   
    \item{\code{Y.nD}}{Person-years among persons without diabetes}
    \item{\code{D.DM}}{Number of deaths among persons with diabetes}	   
    \item{\code{Y.DM}}{Person-years among persons with diabetes}
  }
}
\details{Based on registers of the Danish population. Only included for
  illustrative purposes. Cannot be used as scientifically validaed data.
}
\examples{
data(DMepi)
# Total deaths and person-years in the Danish population
ftable( addmargins( xtabs( cbind( Deaths=D.nD+D.DM,
                                    PYrs=Y.nD+Y.DM ) ~ P + sex,
                           data=DMepi ),
                    2 ),
        row.vars = 1 )
# Deaths and person-years in the population of diabetes patients
round(
ftable( addmargins( xtabs( cbind( Deaths=D.DM,
                                    PYrs=Y.DM ) ~ P + sex,
                           data=DMepi ),
                    2 ),
        row.vars = 1 ) )
}
\keyword{datasets}
