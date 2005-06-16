\name{ewrates}
\alias{ewrates}
\docType{data}
\title{Rates of lung and nasal cancer mortality, and total mortality.}
\description{
England and Wales mortality rates from lung cancer, nasal cancer,
and all causes 1936 - 1980. The 1936 rates are repeated as 1931 rates in
order to accomodate follow up for the \code{\link{nickel}} study.
}
\usage{data(ewrates)}
\format{
  A data frame with 150 observations on the following 5 variables:
\tabular{rl}{
\code{id}: \tab Subject identifier (numeric) \cr
    \code{year}  \tab Calendar period, 1931: 1931--35, 1936: 1936--40,
    \ldots \cr
    \code{age}   \tab Age class: 10: 10--14, 15:15--19, \ldots  \cr
    \code{lung}  \tab Lung cancer mortality rate per 1,000,000 py. \cr
    \code{nasal} \tab Nasal cancer mortality rate per 1,000,000 py. \cr
    \code{other} \tab All cause mortality rate per 1,000,000 py.
}
}
}
\source{
From Breslow and Day, Vol II, Appendix IX.
}
\examples{
data(ewrates)
str(ewrates)
}
\keyword{datasets}
