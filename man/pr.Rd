\name{pr}
\alias{pr}
\docType{data}
\title{Diabetes prevance as of 2010-01-01 in Denmark
}
\description{
Diabetes prevalence as of 2010-01-01 in Denmark in 1-year age classes by sex.
}
\usage{data("pr")}
\format{
  A data frame with 200 observations on the following 4 variables.
  \describe{
    \item{\code{A}}{Numeric, age, 0-99}
    \item{\code{sex}}{Sex, a factor with levels \code{M} \code{F}}
    \item{\code{X}}{Number of diabetes patients}
    \item{\code{N}}{Population size}
  }
}
\examples{
data(pr)
str(pr)
}
\keyword{datasets}
