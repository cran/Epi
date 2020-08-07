\name{births}
\alias{births}
\docType{data}
\title{Births in a London Hospital}
\description{
Data from 500 singleton births in a London Hospital
}
\usage{data(births)}
\format{
  A data frame with 500 observations on the following 8 variables.
  \tabular{rl}{
    \code{id}: \tab Identity number for mother and baby. \cr
    \code{bweight}: \tab Birth weight of baby. \cr
    \code{lowbw}: \tab Indicator for birth weight less than 2500 g. \cr
    \code{gestwks}: \tab Gestation period. \cr
    \code{preterm}: \tab Indicator for gestation period less than 37 weeks. \cr
    \code{matage}: \tab Maternal age. \cr
    \code{hyp}: \tab Indicator for maternal hypertension. \cr
    \code{sex}: \tab Sex of baby: 1:Male, 2:Female. \cr
  }
}
\source{
Anonymous
}
\references{
Michael Hills and Bianca De Stavola (2002). A Short Introduction to
Stata 8 for Biostatistics, Timberlake Consultants Ltd
}
\examples{
data(births)
}
\keyword{datasets}
