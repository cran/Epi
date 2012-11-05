\name{B.dk}
\alias{B.dk}
\docType{data}
\title{Births in Denmark by year and month of birth and sex}
\description{
  The number of live births as entered from printed publications from
  Statistics Denmark.
}
\usage{data(B.dk)}
\format{
  A data frame with 1248 observations on the following 4 variables.
  \describe{
    \item{\code{year}}{Year of birth}
    \item{\code{month}}{Month of birth}
    \item{\code{m}}{Number of male births}
    \item{\code{f}}{Number of female births}
  }
}
\details{
  Division of births by month and sex is only avaialble for the years 1957--69
  and 2002ff. For the remaining period, the total no. births in each month is
  divided between the sexes so that the fraction of boys is equal to the overall
  fraction for the years where the sex information is available.
  
  There is a break in the series at 1920, when Sonderjylland was joined to Denmark.
}
\source{
  Statistiske Undersogelser nr. 19: Befolkningsudvikling og sundhedsforhold
                                    1901-60, Copenhagen 1966.
  Befolkningens bevaegelser 1957.
  Befolkningens bevaegelser 1958.
  ...
  Befolkningens bevaegelser 2003.
  Befolkningens bevaegelser 2004.
  Vital Statistics 2005.
  Vital Statistics 2006.
}
\examples{
data( B.dk )
str( B.dk )
attach( B.dk )
# Plot the no of births and the M/F-ratio
par( las=1, mar=c(4,4,2,4) )
matplot( year+(month-0.5)/12,
         cbind( m, f ),
         bty="n", col=c("blue","red"), lty=1, lwd=1, type="l",
         ylim=c(0,5000),
         xlab="Date of birth", ylab="" )
usr <- par()$usr
mtext( "Monthly no. births in Denmark", side=3, adj=0, at=usr[1], line=1/1.6 )
text( usr[1:2] \%*\% cbind(c(19,1),c(19,1))/20,
      usr[3:4] \%*\% cbind(c(1,19),c(2,18))/20, c("Boys","Girls"), col=c("blue","red"), adj=0 ) 
lines( year+(month-0.5)/12, (m/(m+f)-0.5)*30000, lwd=1 )
axis( side=4, at=(seq(0.505,0.525,0.005)-0.5)*30000, labels=c("","","","",""), tcl=-0.3 )
axis( side=4, at=(50:53/100-0.5)*30000, labels=50:53, tcl=-0.5 )
axis( side=4, at=(0.54-0.5)*30000, labels="\% boys", tick=FALSE, mgp=c(3,0.1,0) )
abline( v=1920, col=gray(0.8) )
}
\keyword{datasets}
