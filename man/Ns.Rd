\name{Ns}
\alias{Ns}
\title{
  Natural splines - (cubic splines linear beyond outermost knots) with
  convenient specification of knots and possibility of centering and
  detrending.
}
\description{
  This function is partly for convenient specification of natural splines
  in practical modeling. The convention used is to take the smallest
  and the largest of the supplied knots as boundary knots. It also has
  the option of centering the effects provided at a chosen reference
  point as well as projecting the columns on the orthogonal space to
  that spanned by the intercept and the linear effect of the variable.
}
\usage{
Ns( x, ref = NULL, df = NULL,
                knots = NULL,
            intercept = FALSE,
       Boundary.knots = NULL,
              detrend = FALSE )
}
\arguments{
  \item{x}{A variable.}
  \item{ref}{Scalar. Reference point on the \code{x}-scale, where the
    resulting effect will be 0.}
  \item{df}{degrees of freedom.}
  \item{knots}{knots to be used both as boundary and internal knots. If
    \code{Boundary.knots} are given, this will be taken as the set of
    internal knots.}
  \item{intercept}{Should the intercept be included in the resulting
    basis? Ignored if any of \code{ref} or \code{detrend} is given.}
  \item{Boundary.knots}{The boundary knots beyond which the spline is
    linear.}
  \item{detrend}{If \code{TRUE}, the columns of the spline basis will be
    projected to the orthogonal of \code{cbind(1,x)}. Optionally
    \code{detrend} can be given as a vector of non-negative numbers used
    to define an inner product as \code{diag(detrend)} for projection on
    the orthogonal to \code{cbind(1,x)}. The default is projection
    w.r.t. the inner product defined by the identity matrix.}
}
\value{
  A matrix of dimension c(length(x),df) where either \code{df} was
  supplied or if \code{knots} were supplied, \code{df = length(knots) -
  intercept}. \code{Ns} returns a spline basis which is centered at
  \code{ref}. \code{Ns} with the argument \code{detrend=TRUE} returns a
  spline basis which is orthogonal to \code{cbind(1,x)} with respect to
  the inner product defined by the positive definite matrix
  \code{diag(weight)} (an assumption which is checked).
}
\author{
  Bendix Carstensen
}
\note{
  The need for this function is primarily from analysis of rates in
  epidemiology and demography, where the dataset are time-split records
  of follow-up, and the range of data therefore rarely is of any
  interest (let alone meaningful). In Poisson modeling of rates
  based on time-split records one should aim at having the same number
  of \emph{events} between knots, rather than the same number of
  observations.
}
\examples{
require(splines)
require(stats)
require(graphics)

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

# Split follow-up in 1-year age intervals
dms <- splitLexis( dml, time.scale="Age", breaks=0:100 )
summary( dms )

# Model  age-specific rates using Ns with 6 knots
# and period-specific RRs around 2000 with 4 knots
# with the same number of deaths between each pair of knots
n.kn <- 6
( a.kn <- with( subset(dms,lex.Xst=="Dead"),
                quantile( Age+lex.dur, probs=(1:n.kn-0.5)/n.kn ) ) )
n.kn <- 4
( p.kn <- with( subset(dms,lex.Xst=="Dead"),
                quantile( Per+lex.dur, probs=(1:n.kn-0.5)/n.kn ) ) )
m1 <- glm( lex.Xst=="Dead" ~ Ns( Age, kn=a.kn ) +
                             Ns( Per, kn=p.kn, ref=2000 ),
           offset = log( lex.dur ), family=poisson, data=dms )

# Plot estimated age-mortality curve for the year 2005 and knots chosen:
nd <- data.frame(Age=40:90,Per=2005,lex.dur=1000)
par( mfrow=c(1,2) )
matplot( nd$Age, ci.pred( m1, newdata=nd ),
         type="l", lwd=c(3,1,1), lty=1, col="black", log="y",
         ylab="Mortality rates per 1000 PY", xlab="Age (years)", las=1 )
rug( a.kn, lwd=2 )
}
\keyword{regression}
