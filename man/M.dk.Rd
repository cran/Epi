\name{M.dk}
\alias{M.dk}
\docType{data}
\title{Mortality in Denmark 1974 ff.}
\description{
  Mortality in one-year classes of age (0-98,99+) and period (1974 ff.) in Denmark.
}
\usage{data(M.dk)}
\format{
  A data frame with 6400 observations on the following 6 variables.
  \describe{
    \item{\code{A}}{Age-class, 0-98, 99:99+}
    \item{\code{sex}}{Sex. 1:males, 2:females}
    \item{\code{P}}{Period (year) of death}
    \item{\code{D}}{Number of deaths}
    \item{\code{Y}}{Number of person-years}
    \item{\code{rate}}{Mortality rate per 1000 person-years}
  }
}
\details{
  Deaths in ages over 100 are in the class labelled 99. Risk time is
  computed by tabulation of the risk time in \code{\link{Y.dk}}, except
  for the class 99+ where the average of the population size in ages
  99+ at the first and last date of the year is used.}
\source{
  \url{http://www.statistikbanken.dk/statbank5a/SelectTable/omrade0.asp?SubjectCode=02&PLanguage=1&ShowNews=OFF}
}
\examples{
data(M.dk)
str(M.dk)

zz <- xtabs( rate ~ sex+A+P, data=M.dk )
zz[zz==0] <- NA # 0s makes log-scale plots crash
par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
for( i in 1:2 )
{
matplot( dimnames(zz)[[2]], zz[i,,],
         lty=1, lwd=1, col=rev(heat.colors(37)),
         log="y", type="l", ylim=range(zz,na.rm=TRUE),
         ylab="", xlab="", yaxt="n" )
text( 0, max(zz,na.rm=TRUE), c("M","F")[i], font=2, adj=0:1, cex=2, col="gray" )
if( i==1 ) axis( side=2, las=1 )
}
mtext( side=1, "Age", line=2, outer=TRUE )
mtext( side=2, "Mortality rate", line=2, outer=TRUE )
}
\keyword{datasets}
