\name{ hivDK }
\alias{ hivDK }
\docType{ data }
\title{ hivDK: seroconversion in a cohort of Danish men}
\description{
  Data from a survey of HIV-positivity of a cohort of Danish
  men followed by regular tests from 1983 to 1989.
  }
\usage{
  data(hivDK)
  }
\format{
  A data frame with 297 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{ID of the person}
    \item{\code{entry}}{Date of entry to the study. Date variable.}
    \item{\code{well}}{Date last seen seronegative. Date variable.}
    \item{\code{ill}}{Date first seen seroconverted. Date variable.}
    \item{\code{bth}}{Year of birth minus 1950.}
    \item{\code{pyr}}{Annual number of sexual partners.}
    \item{\code{us}}{Indicator of wheter the person has visited the USA.}
  }
}
\source{
  Mads Melbye, Statens Seruminstitut.
  }
\references{
  Becker N.G. and Melbye M.: Use of a log-linear model to
  compute the empirical survival curve from interval-censored data,
  with application to data on tests for HIV-positivity, Australian
  Journal of Statistics, 33, 125--133, 1990.

  Melbye M., Biggar R.J., Ebbesen P., Sarngadharan M.G., Weiss
  S.H., Gallo R.C. and Blattner W.A.: Seroepidemiology of HTLV-III
  antibody in Danish homosexual men: prevalence, transmission and
  disease outcome. British Medical Journal, 289, 573--575, 1984.
  }
\examples{
  data(hivDK)
  str(hivDK) 
  }
\keyword{ datasets }
