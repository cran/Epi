\name{simLexis}
\alias{simLexis}
\alias{nState}
\alias{pState}
\alias{plot.pState}
\title{Simulate a Lexis object representing follow-up in a multistate
  model.}

\description{Based on a (pre-)\code{Lexis} object representing persosn
  at given states and tines, and full specification of transition
  intensities between states in the form of fitted Poisson models, this
  function simulates transition times and -types for persons and returns
  a \code{Lexis} object representing the simulted cohort.}

\usage{
simLexis( Tr, init, time.pts = 0:50/2, N = 1,
              lex.id = 1:(N*nrow(init)), type = "glm-mult" )
 nState( obj, at, from, time.scale = 1 )
 pState( nSt, perm = 1:ncol(nSt) )
\method{plot}{pState}( x, col = rainbow(ncol(x)),
                  border = "transparent",
                    xlab = "Time",
                    ylab = "Probability", ... )
   }
\arguments{
  \item{Tr}{A named list of named lists. The names of the lists are
    names of the transient states in the model, and the names of the
    list elements are the names of the states reachable from this. See details.}
  \item{init}{A (pre-)\code{\link{Lexis}} object representing the initial
    state of the persons whose trajectories through the multiple states
    we want to simulate. Must have an attribute "time.since" --- see
    details. Duplicate values of \code{lex.id} is non-sensical and not
    accepted.}
  \item{time.pts}{Numerical vector of times since start. Cumulative
    rates for transitions are computed at these times after stater and
    entry state. Simulation is only done till time \code{max(time.pts)}
    after start, where persons are censored.}
  \item{N}{Numeric. How many persons should be simulated. \code{N}
    persons with covariate configuration of each row in \code{init}
    will be simulated.}
  \item{lex.id}{Vector of ids of the simulated persons. Useful when
    simulating in chunks.}
  \item{type}{Not implemented (yet); \code{simLexis} only works if
    elements of \code{Tr} are glm objects.}
  \item{obj}{A \code{Lexis} object.}
  \item{from}{The point on the time scale \code{time.scale} from which
    we start counting.}
  \item{time.scale}{The timescale to which \code{from} refer.}
  \item{at}{Time points (after \code{from}) where the number of persons
    in each state is to be computed.}
  \item{nSt}{A table obtained by \code{nState}.}
  \item{perm}{A permutation of columns used before cumulating row-wise
    and taking percentages.}
  \item{x}{An object of class \code{pState}, e.g. created by \code{pState}.}
  \item{col}{Colors for filling the areas between curves.}
  \item{border}{Colors for outline of the areas between curves.}
  \item{xlab}{Label on x-axis}
  \item{ylab}{Label on y-axis}
  \item{...}{Further arguments passed on to \code{plot}.}
}
\details{The simulation command \code{simLexis} is not defined as a
  method for \code{Lexis} objects, because the input is not really a
  \code{Lexis} object, the \code{Lexis}-like object is merely
  representing a prevalent population and a specification of which
  variables that are timescales. The variables \code{lex.dur} and
  \code{lex.Xst} are ignored (and overwritten) if present. The core
  input is the list \code{Tr} giving the transitions.

  The components of \code{Tr} represents the transition intensities
  between states. Currently only implemented for \code{type =
  "glm-mult"} which means that they are assumed to be \code{glm}
  objects, specifically Poisson models with \code{log(lex.dur)} as
  offset. Thus the transition from state \code{A} to \code{B}, say, is
  assumed modelled by a glm with Poisson family, log link and an offset
  \code{log(lex.dur)}. The resulting model object is assumed stored in
  \code{Tr$A$B}. Thus names of the elements of \code{Tr} are names of
  transient states, and the names of the elements of these are the names
  of states reachable from these.

  The \code{\link{Lexis}} object \code{init} must contain values of all
  variables used in any of the \code{glm} objects in \code{Tr}. Moreover
  the attributes \code{time.scales} and \code{time.since} must be
  present. The attribute \code{time.since} is a character vector of the
  same length as \code{time.scales} and the elements value \code{"A"} if
  the corresponding time scale is defined as "time since entry into
  state A", otherwise the value is \code{""}. If not present it will be
  set to a vector of \code{""}s, which is only OK if no time scales are
  defined as time since entry to a state.

  The function \code{\link{Lexis}} automatically generates an attribute
  \code{time.since}, and \code{\link{cutLexis}} updates it when new time
  scales are defined. Hence, the simplest way of defining a initial
  (pre-)Lexis object representing a current state of a (set of) persons
  to be followed through a multistate model is to take \code{NULL} rows
  of an existing Lexis object (normally the one used for estimation),
  and so ensuring that all relevant attributes and state levels are
  properly defined. See the example code.

  The prevalence function \code{nState} computes the distribution of
  individuals in different states at prespecified times. Only sensible
  for a simulated \code{Lexis} object. The function \code{pState} takes
  a matrix as output by \code{nState} and computes the row-wise
  cumulative probabilities across states, and leaves an object of class
  \code{pState}, suitable for plotting.
}

\value{\code{simLexis} returns a \code{\link{Lexis}} object representing
  the experience of a population starting as \code{init} followed
  through the states according to the transitions in \code{Tr}.

  The function \code{nState} returns a table of persons classified by
  states at each of the times in \code{at}. Note that this function can
  easily produce meaningless results, for example if applied to a
  \code{Lexis} object not created by simulation. If you apply it to a
  \code{Lexis} object generated by \code{simLexis}, you must make sure
  that you start (\code{from}) the point where you started the
  simulation on the correct timescale. The resulting object has class
  \code{"pState"} and inherits from \code{"matrix"}.
}

\author{Bendix Carstensen, \url{BendixCarstensen.com}.}
\seealso{
  \code{\link{Lexis}},
  \code{\link{cutLexis}},
  \code{\link{splitLexis}}
}
\examples{
data(DMlate)
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit = list(Per=dox),
        exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
               data = DMlate[runif(nrow(DMlate))<0.1,] )
# Split follow-up at insulin, introduce a new timescale,
# and split non-precursor states
dmi <- cutLexis( dml, cut = dml$doins,
                      pre = "DM",
                new.state = "Ins",
                new.scale = "t.Ins",
             split.states = TRUE )
# Split the follow in 1-year intervals for modelling
Si <- splitLexis( dmi, 0:30/2, "DMdur" )
# Define knots
nk <- 4
( ai.kn <- with( subset(Si,lex.Xst=="Ins"),
                 quantile( Age+lex.dur, probs=(1:nk-0.5)/nk ) ) )
( ad.kn <- with( subset(Si,lex.Xst=="Dead"),
                 quantile( Age+lex.dur, probs=(1:nk-0.5)/nk ) ) )
( di.kn <- with( subset(Si,lex.Xst=="Ins"),
                 quantile( DMdur+lex.dur, probs=(1:nk-0.5)/nk ) ) )
( dd.kn <- with( subset(Si,lex.Xst=="Dead"),
                 quantile( DMdur+lex.dur, probs=(1:nk-0.5)/nk ) ) )
( td.kn <- with( subset(Si,lex.Xst=="Dead(Ins)"),
                 quantile( t.Ins+lex.dur, probs=(1:nk-0.5)/nk ) ) )

# Fit Poisson models to transition rates
library( splines )
DM.Ins <- glm( (lex.Xst=="Ins") ~ ns( Age  , knots=ai.kn[2:(nk-1)], Bo=ai.kn[c(1,nk)] ) +
                                  ns( DMdur, knots=di.kn[2:(nk-1)], Bo=di.kn[c(1,nk)] ) +
                                  I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
DM.Dead <- glm( (lex.Xst=="Dead") ~ ns( Age  , knots=ad.kn[2:(nk-1)], Bo=ad.kn[c(1,nk)] ) +
                                    ns( DMdur, knots=dd.kn[2:(nk-1)], Bo=dd.kn[c(1,nk)] ) +
                                    I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
Ins.Dead <- glm( (lex.Xst=="Dead(Ins)") ~ ns( Age  , knots=ad.kn[2:(nk-1)], Bo=ad.kn[c(1,nk)] ) +
                                          ns( DMdur, knots=dd.kn[2:(nk-1)], Bo=dd.kn[c(1,nk)] ) +
                                          ns( t.Ins, knots=td.kn[2:(nk-1)], Bo=td.kn[c(1,nk)] ) +
                                          I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="Ins") )

# Stuff the models into an object representing the transitions
Tr <- list( "DM" = list( "Ins"       = DM.Ins,
                         "Dead"      = DM.Dead  ),
           "Ins" = list( "Dead(Ins)" = Ins.Dead ) )
lapply( Tr, names )

# Define an initial object - note the combination of "select=" and NULL
# which ensures that the relevant attributes from the Lexis object 'Si'
# are carried over to 'ini':
ini <- subset(Si,select=1:9)[NULL,]
ini[1:2,"lex.Cst"] <- "DM"
ini[1:2,"Per"] <- 1995
ini[1:2,"Age"] <- 60
ini[1:2,"DMdur"] <- 5
ini[1:2,"sex"] <- c("M","F")
str(ini)

# Simulate 200 of each sex using the estimated model
simL <- simLexis( Tr, ini, time.pts=seq(0,50,0.5), N=200 )
summary( simL, by="sex" )

# Give the number of persons in each state at a set of times
nSt <- nState( subset(simL,sex=="M"),
               at=seq(0,15,0.2), from=1995, time.scale="Per" )
nSt

# Show the cumulative prevalences in a different order than that of the
# state-level ordering and plot them using all defaults
pp <- pState( nSt, perm=c(1,2,4,3) )
head( pp )
plot( pp )

# A more useful set-up of the graph
clr <- c("orange2","forestgreen")
par( las=1 )
plot( pp, col=clr[c(2,1,1,2)] )
lines( as.numeric(rownames(pp)), pp[,2], lwd=2 )
mtext( "60 year old male, diagnosed 1995", side=3, line=2.5, adj=0 )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[1] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[2] )
axis( side=4 )
}
\keyword{survival}
