\name{testisDK}
\alias{testisDK}
\docType{data}
\title{Testis cancer incidence in Denmark, 1943--1996
}
\description{Number of testiscancer cases and male person-years in the
  Danish population 1943--1996}
\usage{data(testisDK)}
\format{
  A data frame with 4860 observations on the following 4 variables.
  \describe{
    \item{\code{A}}{Age class, 0,1,2,...,89}
    \item{\code{P}}{Year, 1943,...,1996}
    \item{\code{D}}{Number of testis cancer cases}
    \item{\code{Y}}{Person years}
  }
}
\source{The Danish Cancer Registry
}
\examples{
data(testisDK)
head(testisDK)
}
\keyword{datasets}
