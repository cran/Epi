\name{Ns}
\alias{Ns}
\title{Natural splines - cubic splines, made linear beyond outermost knots.}
\description{This function is just for convenient specification of
  natural splines in practical modelling. The convention used is to take
  the smallest and the largest of the supplied knots as boundary knots.}
\usage{
Ns(x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = NULL)
}
\arguments{
  \item{x}{A variable.}
  \item{df}{degrees of freedom.}
  \item{knots}{knots to be used both as boundary and internal knots. If
    \code{Boundary.knots} are given, this will be taken as the set of
    internal knots.}
  \item{intercept}{Should the intercept be included in the resulting basis?}
  \item{Boundary.knots}{The boundary knots beyond which the spline is linear.}
}
\details{The function is a convenience wrapper, which gives the same
  result as \code{ns}, but allows omission of \code{Boundary.knots}. It
  has the same arguments as \code{\link[splines:splines]{ns}}, and in
  the same order.}
\value{ A matrix of dimension c(length(x),df) where either \code{df} was
  supplied or if \code{knots} were supplied,
  \code{df = length(knots) - intercept}. }
\author{Bendix Carstensen}
\note{The need for this function is from analysis of rates in
  epidemiology and demography, where the dataset are time-split records
  of follow-up, and the range of data therefore rarely is of any
  interest (let alone is meaningful). In Poisson modelling of rates based
  on time-split records one should aim at having the same number of
  \emph{events} between knots, rather than the same number of
  observations.}
\examples{
require(splines); require(stats); require(graphics)
ns( women$height, df = 3)
Ns( women$height, knots=c(63,59,71,67) )
# Gives the same results as ns:
summary( lm(weight ~ ns(height, df = 3), data = women) )
summary( lm(weight ~ Ns(height, df = 3), data = women) )

# Get the diabetes data and set up as Lexis object
data(DMlate)
DMlate <- DMlate[sample(1:nrow(DMlate),500),]
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit = list(Per=dox),
        exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
               data = DMlate )

# Split in 1-year age intervals
dms <- splitLexis( dml, time.scale="Age", breaks=0:100 )
summary( dms )

# Model age-specific rates using Ns with 7 knots
n.kn <- 7
( a.kn <- with( subset(dms,lex.Xst=="Dead"),
                quantile( Age+lex.dur, probs=(1:n.kn-0.5)/n.kn ) ) )
m0 <- glm( lex.Xst=="Dead" ~ Ns( Age, kn=a.kn ),
           offset = log( lex.dur ), family=poisson, data=dms )

# Plot estimated curve:
matplot( 50:100, ci.exp( m0, ctr.mat=cbind(1,Ns(50:100,knots=a.kn)) ),
         type="l", lwd=c(3,1,1), lty=1, col="black", log="y",
         ylab="Mortality rates per 1000 PY", xlab="Age (years)", las=1 )
rug( a.kn, lwd=2 )
}
\keyword{regression}
