\name{N.dk}
\alias{N.dk}
\docType{data}
\title{Population size in Denmark}
\description{
  The population size at 1st January in ages 0-99.
}
\usage{data(N.dk)}
\format{
  A data frame with 7200 observations on the following 4 variables.
  \describe{
    \item{\code{sex}}{Sex, 1:males, 2:females}
    \item{\code{A}}{Age. 0:0, 1:1, ..., 98:98, 99:99+}
    \item{\code{P}}{Year}
    \item{\code{N}}{Number of persons alive at 1st January year \code{P}}
  }
}
\source{
  \url{http://www.statistikbanken.dk/statbank5a/SelectTable/omrade0.asp?SubjectCode=02&PLanguage=1&ShowNews=OFF}
}
\examples{
data(N.dk)
str(N.dk)
with(N.dk,addmargins(tapply(N,list(P,sex),sum),2))
with(subset(N.dk,P==max(P)),addmargins(tapply(N,list(A,sex),sum)))
}
\keyword{datasets}
