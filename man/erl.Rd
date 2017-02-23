\name{erl}
\alias{surv1}
\alias{surv2}
\alias{erl1}
\alias{erl}
\alias{yll}
\title{Compute survival functions from rates and expected residual
  lifetime in an illness-death model as well as years of life lost to disease.
}
\description{
  These functions compute survival functions from a set of mortality and
  disease incidence rates in an illness-death model. Expected residual
  life time can be computed under various scenarios by the \code{erl}
  function, and areas between survival functions can be computed under
  various scenarios by the \code{yll} function. Rates are assumed
  supplied for equidistant intervals of length \code{int}.
  }
\usage{
  surv1( int, mu ,                age.in = 0, A = NULL )
   erl1( int, mu ,                age.in = 0 ) 
  surv2( int, muW, muD, lam,      age.in = 0, A = NULL )
    erl( int, muW, muD, lam=NULL, age.in = 0, A = NULL,
         immune = is.null(lam), yll=TRUE, note=TRUE )
    yll( int, muW, muD, lam=NULL, age.in = 0, A = NULL,
         immune = is.null(lam), note=TRUE )
  }
\arguments{
  \item{int}{
  Scalar. Length of intervals that rates refer to.    
  }
  \item{mu}{
  Numeric vector of mortality rates at midpoints of intervals of length \code{int}
  }
  \item{muW}{
  Numeric vector of mortality rates among persons in the "Well" state  at
  midpoints of intervals of length \code{int}. Left endpoint of first
  interval is \code{age.in}. 
  }
  \item{muD}{
  Numeric vector of mortality rates among persons in the "Diseased" state
  at midpoints of intervals of length \code{int}.  Left endpoint of first
  interval is \code{age.in}. 
  }
  \item{lam}{
  Numeric vector of disease incidence rates among persons in the "Well" state
  at midpoints of intervals of length \code{int}. Left endpoint of first
  interval is \code{age.in}. 
  }
  \item{age.in}{
  Scalar indicating the age at the left endpoint of the first interval. 
  }
  \item{A}{
  Numeric vector of conditioning ages for calculation of survival
  functions.   
  }
  \item{immune}{
  Logical. Should the years of life lost to the disease be computed
  using assumptions that non-diseased individuals are immune to the
  disease (\code{lam}=0) and that their mortality is yet still \code{muW}.
  }
  \item{note}{
  Logical. Should a warning of silly assumptions be printed?
  }
  \item{yll}{
  Logical. Should years of life lost be included in the result?
  }
}
\details{
  The mortality rates given are supposed to refer to the ages
  \code{age.in+(i-1/2)*int}, \code{i=1,2,3,...}.
  
  The units in which \code{int} is given must correspond to the units in
  which the rates \code{mu}, \code{muW}, \code{muD} and \code{lam} are
  given. Thus if \code{int} is given in years, the rates must be given
  in the unit of events per year.

  The ages in which the survival curves are computed are from
  \code{age.in} and then at the end of \code{length(muW)}
  (\code{length(mu)}) intervals each of length \code{int}.
  
  The \code{age.in} argument is merely a device to account for rates
  only available from a given age. It has two effects, one is that
  labeling of the interval endpoint is offset by this quantity, thus
  starting at \code{age.in}, and the other that the conditioning ages
  given in the argument \code{A} will refer to the ages defined by this.
  
  The \code{immune} argument is \code{FALSE} whenever the disease
  incidence rates are supplied. If set to \code{TRUE}, the years of life
  lost is computed under the assumption that individuals without the
  disease at a given age are immune to the disease in the sense that the
  disease incidence rate is 0, so transitions to the diseased state
  (with presumably higher mortality rates) are assumed not to
  occur. This is a slightly peculiar assumption (but presumably the most
  used in the epidemiological literature) and the resulting object is
  therefore given an attribute, \code{NOTE}, that point this out. The
  default of the \code{surv2} function is to take the possibility of
  disease into account in order to potentially rectify this.}

\value{\code{surv1} and \code{surv2} return a matrix whose first column
  is the ages at the ends of the 
  intervals, thus with \code{length(mu)+1} rows. The following columns
  are the survival functions (since \code{age.in}), and conditional on
  survival till ages as indicated in \code{A}, thus a matrix with
  \code{length(A)+2} columns. Columns are labeled with the actual
  conditioning ages; if \code{A} contains values that are not among the
  endpoints of the intervals used, the nearest smaller interval border
  is used as conditioning age, and columns are named accordingly.

  \code{surv1} returns the survival function for a simple model with one
  type of death, occurring at intensity \code{mu}.
  
  \code{surv2} returns the survival function for a person in the "Well"
  state of an illness-death model, taking into account that the person
  may move to the "Diseased" state, thus requiring all three transition
  rates to be specified. The conditional survival functions are
  conditional on being in the "Well" state at ages given in \code{A}.

  \code{erl1} returns a three column matrix with columns \code{age},
  \code{surv} (survival function) and \code{erl} (expected residual life
  time) with \code{length(mu)+1} rows. 
  
  \code{erl} returns a two column matrix, columns labeled "Well" and
  "Dis", and with row-labels \code{A}. The entries are the expected
  residual life times given survival to \code{A}. If \code{yll=TRUE} the
  difference between the columns is added as a
  third column, labeled "YLL".
  }
\author{Bendix Carstensen, \email{b@bxc.dk}
  }
\seealso{
\code{\link{ci.cum}}
}
\examples{
library( Epi )
data( DMlate )
# Naive Lexis object
Lx <- Lexis( entry = list( age = dodm-dobth ),
              exit = list( age = dox -dobth ),
       exit.status = factor( !is.na(dodth), labels=c("DM","Dead") ),
              data = DMlate )
# Cut follow-up at insulin inception
Lc <- cutLexis( Lx, cut = Lx$doins-Lx$dob,
              new.state = "DM/ins",
       precursor.states = "DM" )
summary( Lc )
# Split in small age intervals
Sc <- splitLexis( Lc, breaks=seq(0,120,1) )
summary( Sc )

# Overview of object
boxes( Sc, boxpos=TRUE, show.BE=TRUE, scale.R=100 )

# Knots for splines
a.kn <- 2:9*10

# Mortality among DM
mW <- glm( lex.Xst=="Dead" ~ Ns( age, knots=a.kn ),
           offset = log(lex.dur),
           family = poisson,
             data = subset(Sc,lex.Cst=="DM") )

# Mortality among insulin treated
mI <- update( mW, data = subset(Sc,lex.Cst=="DM/ins") )

# Total motality
mT <- update( mW, data = Sc )

# Incidence of insulin inception
lI <- update( mW, lex.Xst=="DM/ins" ~ . )

# From these we can now derive the fitted rates in intervals of 1 year's
# length. In real applications you would use much smaller interval like
# 1 month:
# int <- 1/12 
int <- 1

# Prediction frame to return rates in units of cases per 1 year
# - we start at age 40 since rates of insulin inception are largely
# indeterminate before age 40
nd <- data.frame( age = seq( 40+int, 110, int ) - int/2,
              lex.dur = 1 )
muW <- predict( mW, newdata = nd, type = "response" )
muD <- predict( mI, newdata = nd, type = "response" )
lam <- predict( lI, newdata = nd, type = "response" )

# Compute the survival function, and the conditional from ages 50 resp. 70
s1 <- surv1( int, muD, age.in=40, A=c(50,70) )
round( s1, 3 )

s2 <- surv2( int, muW, muD, lam, age.in=40, A=c(50,70) )
round( s2, 3 )

# How much is YLL overrated by ignoring insulin incidence?
round( YLL <- cbind(
yll( int, muW, muD, lam, A = 41:90, age.in = 40 ),
yll( int, muW, muD, lam, A = 41:90, age.in = 40, immune=TRUE ) ), 2 )[seq(1,51,10),]

par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6, bty="n", las=1 )
matplot( 40:90, YLL,
         type="l", lty=1, lwd=3,
         ylim=c(0,10), yaxs="i", xlab="Age" )
}
\keyword{survival}
