\name{DMconv}
\alias{DMconv}
\docType{data}
\title{Conversion to diabetes}
\description{
  Data from a randomized intervention study ("Addition") where persons with
  prediabetic conditions are followed up for conversion to diabetes (DM).
  Conversion dates are interval censored.
  Original data are not published yet, so id-numbers have been changed and
  all dates have been randomly perturbed.
  }
\usage{data(DMconv)}
 \format{
  A data frame with 1519 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{Person identifier}
    \item{\code{doe}}{Date of entry, i.e. first visit.}
    \item{\code{dlw}}{Date last seen well, i.e. last visit without DM.}
    \item{\code{dfi}}{Date first seen ill, i.e. first visit with DM.}
    \item{\code{gtol}}{Glucose tolerance. Factor with levels:
        1="IFG" (impaired fasting glucose), 2="IGT" (impaired glucose tolerance).}
    \item{\code{grp}}{Randomization. Factor with levels:
        1="Intervention", 2="Control".}
  }
}
\source{
  Signe Saetre Rasmussen, Steno Diabetes Center. The Addition Study.
}
\examples{
data(DMconv)
str(DMconv)
head(DMconv)
}
\keyword{datasets}
