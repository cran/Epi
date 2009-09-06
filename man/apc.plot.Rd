\name{apc.plot}
\alias{apc.plot}
\title{Plot the estimates from a fitted Age-Period-Cohort model}
\description{
This function plots the estimates created by \code{\link{apc.fit}} in a single 
graph. It just calls \code{\link{apc.frame}} after computing some sensible 
values of the parameters, and subsequently plots the estimates using 
\code{\link{apc.lines}}.
}
\usage{
apc.plot(obj, r.txt = "Rate", ...)
}
\arguments{
  \item{obj}{ An object of class \code{apc}. }
  \item{r.txt}{ The text to put on the vertical rate axis. }
  \item{\dots}{ Additional arguments passed on to \code{\link{apc.lines}}. }
}
\value{ 
 A numerical vector of length two, with names
  \code{c("cp.offset","RR.fac")}. The first is the offset for the cohort
  period-axis, the second the multiplication factor for the rate-ratio
  scale. Therefore, if you want to plot at \code{(x,y)} in the right panel,
  use \code{(x-res["cp.offset"],y/res["RR.fac"])=(x-res[1],y/res[2])}.
  This vector should be supplied for the parameter \code{frame.par} to
  \code{\link{apc.lines}} if more sets of estimates is plotted in the
  same graph.
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://www.pubhealth.ku.dk/~bxc}
} 
\seealso{
  \code{\link{apc.lines}} ,\code{\link{apc.frame}}, \code{\link{apc.fit}}
}
\examples{
data( lungDK )
attach( lungDK )
apc1 <- apc.fit( A=Ax, P=Px, D=D, Y=Y/10^5 )
fp <- apc.plot( apc1 )
apc.lines( apc1, frame.par=fp, drift=1.01, col="red" )
for( i in 1:11 )
  apc.lines( apc1, frame.par=fp, drift=1+(i-6)/100, col=rainbow(12)[i] )
}
\keyword{hplot}
