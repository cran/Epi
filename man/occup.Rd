\name{occup}
\alias{occup}
\docType{data}
\title{ A small occupational cohort }
\description{This is the data that is behind the illustrative Lexis
diagram in Breslow & Day's book on case-control studies.}
\usage{data(occup)}
\format{
  A data frame with 13 observations on the following 4 variables.
  \describe{
    \item{\code{AoE}}{a numeric vector, Age at Entry}
    \item{\code{DoE}}{a numeric vector, Date of entry}
    \item{\code{DoX}}{a numeric vector, Date of eXit}
    \item{\code{Xst}}{eXit status \code{D}-event, \code{W}-withdrawal,
                     \code{X}-censoring}
  }
}
\references{
Breslow & Day: Statistical Methods in Cancer Research, vol 1:
The analysis of case-control studies, figure 2.2, p. 48.}
\examples{
data(occup)
lx <- Lexis( entry = list( per=DoE, age=AoE ),
              exit = list( per=DoX ),
      entry.status = "W",
       exit.status = Xst,
              data = occup )
plot( lx )
# Split follow-up in 5-year classes
sx <- splitLexis( lx, seq(1940,1960,5), "per" )
sx <- splitLexis( sx, seq(  40,  60,5), "age" )
plot( sx )

# Plot with a bit more paraphernalia and a device to get
# the years on the same physical scale on both axes
ypi <- 2.5 # Years per inch
x11( height=15/ypi+1, width=20/ypi+1 )    # add an inch in each direction for
par( mai=c(3,3,1,1)/4, mgp=c(3,1,0)/1.6 ) # the margins set in inches by mai=
plot(sx,las=1,col="black",lty.grid=1,lwd=2,type="l",
     xlim=c(1940,1960),ylim=c(40,55),xaxs="i",yaxs="i",yaxt="n",
     xlab="Calendar year", ylab="Age (years)")
axis( side=2, at=seq(40,55,5), las=1 )
points(sx,pch=c(NA,16)[(sx$lex.Xst=="D")+1] )
box()
# Annotation with the person-years
PY.ann.Lexis( sx, cex=0.8 )
}
\keyword{datasets}
