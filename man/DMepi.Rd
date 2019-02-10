\name{DMepi}
\alias{DMepi}
\docType{data}
\title{Epidemiological rates for diabetes in Denmark 1996--2015}
\description{Register based counts and person-years for incidence of
  diabetes and mortality with and without diabetes.
}
\usage{data("DMepi")}
\format{
  A data frame with 4200 observations on the following 8 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{M}, \code{F}}
    \item{\code{A}}{Age class, 0--99}
    \item{\code{P}}{Calendar year, 1996--2016}
    \item{\code{X}}{Number of new diagnoses of diabetes among persons without diabetes}
    \item{\code{D.nD}}{Number of deaths among persons without diabetes}	   
    \item{\code{Y.nD}}{Person-years among persons without diabetes}
    \item{\code{D.DM}}{Number of deaths among persons with diabetes}	   
    \item{\code{Y.DM}}{Person-years among persons with diabetes}
  }
}
\details{Based on registers of the Danish population. Only included for
  illustrative purposes. Cannot be used as scientifically validated
  data, since small numbers are randomly permuted between units.
}
\examples{
data(DMepi)
# Total deaths and person-years in the Danish population
ftable( addmargins( xtabs( cbind( Deaths=D.nD+D.DM,
                                    PYrs=Y.nD+Y.DM ) ~ P + sex,
                           data=DMepi ),
                    2 ),
        row.vars = 1 )
# Deaths and person-years in the population of diabetes patients
round(
ftable( addmargins( xtabs( cbind( Deaths=D.DM,
                                    PYrs=Y.DM ) ~ P + sex,
                           data=DMepi ),
                    2 ),
        row.vars = 1 ) )

# Model for age-specific incidence rates;
inc <- glm( X ~ sex + Ns( A, knots=seq(30,80,10) ) + P,
                offset = log(Y.nD),
                family = poisson,
                  data = DMepi )

# Predict for men and women separately in 2010:
ndm <- data.frame( sex="M", A=20:90, P=2010, Y.nD=1000 )
ndf <- data.frame( sex="F", A=20:90, P=2010, Y.nD=1000 )
prM <- ci.pred( inc, ndm )
prF <- ci.pred( inc, ndf )
matplot( ndm$A, cbind(prM,prF),
         type="l", lty=1, lwd=c(3,1,1),
         col=rep(c("blue","red"),each=3),
         log="y", xlab="Age", ylab="DM incidence per 1000 PY" )

# This is a proportional hazards model - add sex-age interaction
int <- update( inc, . ~ . + sex:Ns( A, knots=seq(30,80,10) ) )
prM <- ci.pred( int, ndm )
prF <- ci.pred( int, ndf )
matplot( ndm$A, cbind(prM,prF),
         type="l", lty=1, lwd=c(3,1,1),
         col=rep(c("blue","red"),each=3),
         log="y", xlab="Age", ylab="DM incidence per 1000 PY" )

# The rate-ratio is teased out using the ci.exp:
RRp <- ci.exp( inc, list(ndm,ndf) )
RRi <- ci.exp( int, list(ndm,ndf) )

# and added to the plot
matlines( ndm$A, cbind(RRi,RRp),
          type="l", lty=1, lwd=c(3,1,1), col=gray(rep(c(0.3,0.7),each=3)) )
abline(h=1)
axis(side=4)
mtext( "Male/Female IRR", side=4, line=2 )
}
\keyword{datasets}
