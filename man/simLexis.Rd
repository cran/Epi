\name{simLexis}
\alias{simLexis}
\alias{nState}
\alias{pState}
\alias{plot.pState}
\alias{lines.pState}
\title{Simulate a Lexis object representing follow-up in a multistate
  model.}

\description{Based on a (pre-)\code{Lexis} object representing persons
  at given states and times, and full specification of transition
  intensities between states in the form of models for the transition
  rates, this function simulates transition times and -types for persons
  and returns a \code{Lexis} object representing the simulated
  cohort. The simulation scheme accommodates multiple timescales,
  including time since entry into an intermediate state, and accepts
  fitted Poisson models, Cox-models or just a function as specification
  of rates.}

\usage{
simLexis( Tr, init,
           N = 1,
      lex.id,
     t.range = 20,
       n.int = 101,
    time.pts = seq(0,t.range,length.out=n.int) )
nState( obj, at, from, time.scale = 1 )
pState( nSt, perm = 1:ncol(nSt) )
\method{plot}{pState}( x,
                     col = rainbow(ncol(x)),
                  border = "transparent",
                    xlab = "Time",
                    ylim = 0:1,
                    ylab = "Probability", ... )
\method{lines}{pState}( x,
                      col = rainbow(ncol(x)),
                   border = "transparent", ... )
   }
\arguments{
  \item{Tr}{A named list of named lists. The names of the list are names
    of the transient states in the model. Each list element is again a
    named list. The names of the elements of this inner list are the
    names of the states reachable from the state with name equal to the
    list. Elements of the intter lists represent transitions. See details.}
  \item{init}{A (pre-)\code{\link{Lexis}} object representing the initial
    state of the persons whose trajectories through the multiple states
    we want to simulate. Must have attributes "time.scales" and "time.since" --- see
    details. Duplicate values of \code{lex.id} are not sensible and not
    accepted.}
  \item{N}{Numeric. How many persons should be simulated. \code{N}
    persons with covariate configuration of each row in \code{init}
    will be simulated. Either a scalar or a vector of length \code{nrow(init)}.}
  \item{lex.id}{Vector of ids of the simulated persons. Useful when
    simulating in chunks.}
  \item{t.range}{Numerical scalar. The range of time over which to
    compute the cumulative rates when simulating. Simulted times
    beyond this will result in an obervation censored at \code{t.range}
    after entry.}
  \item{n.int}{Number of intervals to use when computing (cumulative) rates.}
  \item{time.pts}{Numerical vector of times since start. Cumulative
    rates for transitions are computed at these times after stater and
    entry state. Simulation is only done till time \code{max(time.pts)}
    after start, where persons are censored. Must start with 0.}
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
  \item{ylim}{Limits on y-axis}
  \item{ylab}{Label on y-axis}
  \item{...}{Further arguments passed on to \code{plot}.}
}
\details{The simulation command \code{simLexis} is not defined as a
  method for \code{Lexis} objects, because the input is not a
  \code{Lexis} object, the \code{Lexis}-like object is merely
  representing a prevalent population and a specification of which
  variables that are timescales. The variables \code{lex.dur} and
  \code{lex.Xst} are ignored (and overwritten) if present. The core
  input is the list \code{Tr} giving the transitions.

  The components of \code{Tr} represents the transition intensities
  between states. The transition from state \code{A} to \code{B}, say,
  is assumed stored in \code{Tr$A$B}. Thus names of the elements of
  \code{Tr} are names of transient states, and the names of the elements
  of each these are the names of states reachable from the corresponding
  transient state.

  The transition intensities are assumed modelled by either a glm with
  Poisson family or a Cox-model. In both cases the timescale(s) in the
  model must be using the names fo the timescales in a Lexis object
  representng the follow-up in a cohort, and the risk time must be taken
  from the variable \code{lex.dur} --- see the example.

  Alternatively, an element in \code{Tr} could be a function
  that from a data frame produces transition rates, or specifically
  cumulative transition rates over intervals of length \code{lex.dur}.

  The pre-\code{Lexis} object \code{init} must contain values of all
  variables used in any of the objects in \code{Tr}, as well as all
  timescales - even those not used in the models. Moreover, the
  attributes \code{time.scales} and \code{time.since} must be
  present. The attribute \code{time.since} is a character vector of the
  same length as \code{time.scales} and an element has value \code{"A"}
  if the corresponding time scale is defined as
  "time since entry into state \code{A}", otherwise the value is
  \code{""}. If not present it will be set to a vector of \code{""}s,
  which is only OK if no time scales are defined as time since entry to
  a state.

  Note that the variables pre-\code{Lexis} object \code{init} must have
  the same mode and class as in the dataset used for fitting the models
  --- hence the indexing of rows by brackets in the assignment of values used in
  the example below - this way the variables have their attributes
  preserved; using \code{init[,"var"] <-} or \code{init$var <-} replaces
  the variable, whereas  \code{init[1:4,"var"] <-} or
  \code{init$var[1:4] <-} replaces values only and prevents you from
  entering non-existing factor levels etc.
  
  The function \code{\link{Lexis}} automatically generates an attribute
  \code{time.since}, and \code{\link{cutLexis}} updates it when new time
  scales are defined. Hence, the simplest way of defining a initial
  pre-\code{Lexis} object representing a current state of a (set of) persons
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
  simulation on the correct timescale, and you will get funny results if
  you try to tabulate beyond the censoring time for the simulation.
  The resulting object has class \code{"table"}.

  The result from using \code{pState} on the result from \code{nState}
  has class \code{c("pState","matrix")}.
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
DM.Ins <- glm( (lex.Xst=="Ins") ~ Ns( Age  , knots=ai.kn ) +
                                  Ns( DMdur, knots=di.kn ) +
                                  I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
DM.Dead <- glm( (lex.Xst=="Dead") ~ Ns( Age  , knots=ad.kn ) +
                                    Ns( DMdur, knots=dd.kn ) +
                                    I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
Ins.Dead <- glm( (lex.Xst=="Dead(Ins)") ~ Ns( Age  , knots=ad.kn ) +
                                          Ns( DMdur, knots=dd.kn ) +
                                          Ns( t.Ins, knots=td.kn ) +
                                          I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="Ins") )

# Stuff the models into an object representing the transitions
Tr <- list( "DM" = list( "Ins"       = DM.Ins,
                         "Dead"      = DM.Dead  ),
           "Ins" = list( "Dead(Ins)" = Ins.Dead ) )
lapply( Tr, names )

# Define an initial object - note the subsetting that ensures that
# all attributes are carried over
ini <- Si[1,1:9][-1,]
ini[1:2,"lex.Cst"] <- "DM"
ini[1:2,"Per"] <- 1995
ini[1:2,"Age"] <- 60
ini[1:2,"DMdur"] <- 5
ini[1:2,"sex"] <- c("M","F")
str(ini)

# Simulate 200 of each sex using the estimated models in Tr
simL <- simLexis( Tr, ini, time.pts=seq(0,11,0.5), N=200 )
summary( simL )

# Find the number of persons in each state at a set of times.
# Note that the times are shirter than the time-span simulated.
nSt <- nState( subset(simL,sex=="M"),
               at=seq(0,10,0.1), from=1995, time.scale="Per" )
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

# Using a Cox-model for the mortality rates assuming the two mortality
# rates to be proportional:
# When we fit a Cox-model, lex.dur must be used in the Surv() function,
# and the I() constrction must be used when specifying intermediate
# states as covariates, since factors with levels not present in the
# data will create NAs in the parameter vector returned by coxph, which
# in return will crash the simulation machinery.
library( survival )
Cox.Dead <- coxph( Surv( DMdur, DMdur+lex.dur,
                         lex.Xst \%in\% c("Dead(Ins)","Dead")) ~
                   Ns( Age-DMdur, knots=ad.kn ) +
                   I(lex.Cst=="Ins") +
                   I(Per-2000) + sex,
               data = Si )
Cr <- list( "DM" = list( "Ins"       = DM.Ins,
                         "Dead"      = Cox.Dead  ),
           "Ins" = list( "Dead(Ins)" = Cox.Dead ) )
simL <- simLexis( Cr, ini, time.pts=seq(0,11,0.2), N=200 )
summary( simL )
nSt <- nState( subset(simL,sex=="M"),
               at=seq(0,10,0.2), from=1995, time.scale="Per" )
pp <- pState( nSt, perm=c(1,2,4,3) )
plot( pp )
}
\keyword{survival}
