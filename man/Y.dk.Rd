\name{Y.dk}
\alias{Y.dk}
\docType{data}
\title{Population risk time in Denmark}
\description{
  Risk time (person-years) in the Danish population, classified by sex, age, period and
  date of birth in 1-year classes. This corresponds to triangles in a Lexis diagram.
}
\usage{data(Y.dk)}
\format{
  A data frame with 13860 observations on the following 6 variables.
  \describe{
    \item{\code{sex}}{Sex. 1:males, 2:females}
    \item{\code{A}}{One-year age class}
    \item{\code{P}}{Period}
    \item{\code{C}}{Birth cohort}
    \item{\code{Y}}{Person-years}
    \item{\code{upper}}{Indicator of upper triangle in the Lexis diagram}
  }
}
\details{
  The risk time is computed from the population size figures in \code{\link{N.dk}},
  using the formulae devised in B. Carstensen:
  "Demography and epidemiology: Age-Period-Cohort models in the computer age",
  \url{http://www.pubhealth.ku.dk/bs/publikationer/rr-06-1.pdf/}
}
\source{
  \url{http://www.statistikbanken.dk/statbank5a/SelectTable/omrade0.asp?SubjectCode=02&PLanguage=1&ShowNews=OFF}
}
\examples{
data(Y.dk)
str(Y.dk)
# Compute mean age, period for the triangles
attach( Y.dk )
age <- A + (1+upper)/3
per <- P + (2-upper)/3
# Plot a Lexis diagram
library( Epi )
Lexis.diagram( age=c(0,10), date=c(1990,2000), coh.grid=TRUE, int=1 )
box()
# Print the person-years for males there
text( per[sex==1], age[sex==1],
      formatC( Y[sex==1]/1000, format="f", digits=1 ) )
}
\keyword{datasets}
