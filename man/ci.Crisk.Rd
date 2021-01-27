\name{ci.Crisk}
\alias{ci.Crisk}
\alias{sim2ci.Crisk}
\alias{sim2ci.Srisk}
\alias{sim2ci.Stime}
\title{
Compute cumulative risk and/or expected sojourn time from competing
risks models. 
}
\description{
Consider a list of parametric models for rates of competing events, such
as different causes of death, A, B, C, say. From estimates of the
cause-specific rates we can then by simple numerical integration compute
the cumulative risk of being in each state ('Surv' (=no event) and A, B
and C) at different times, as well as the stacked cumulative rates such
as A, A+C, A+C+Surv. Finally, we can compute the expected (truncated)
sojourn times in each state up to each time point.

This function does this for simulated samples from the parameter
vectors of supplied model objects, and computes  the mentioend
quantities with simulation-based confidence intervals. Some call this a
prametric bootstrap.

The times and other covariates determining the cause-specific rates must
be supplied in a data frame which will be used for predicting rates for
all transitions.
}
\usage{
ci.Crisk(mods,
           nd,
          int = mean(diff(nd[, 1])),
           nB = 1000,
         perm = length(mods):0 + 1,
        alpha = 0.05, 
       sim.res = 'none')
sim2ci.Crisk(probs, alpha = 0.05)
sim2ci.Srisk(probs, perm = 1:dim(probs)[2],
                    alpha = 0.05)
sim2ci.Stime(probs, int = attr(probs, "int"),
                    alpha = 0.05)
}
\arguments{
  \item{mods}{A named list of \code{glm}/\code{gam} model objects
    representing the cause-specific rates. If the list is not named the
    function will crash. The names will be used as names for the states
    (competing risks), while the state without any event will be called
    "\code{Surv}".
}  
  \item{nd}{A data frame of prediction points and covariates. Must
    represent midpoints of equidistant intervals.
}
  \item{int}{Numeric, the length of the intervals. Defaults to the
    differences in the first column of \code{nd}.
}
  \item{nB}{Scalar. The number of simulations, that is samples from
    the (posterior) distribution of the model parameters. 
}
\item{perm}{Numerical vector of length \code{lengh(mods)+1} indicating
    the order in which states are to be stacked. The \code{'Surv'} state
    is taken to be the first, the remaining in the reverse order supplied
    in the \code{mods} argument. The default is therefore to stack with
    the survival as the first, which may not be what you normally want.   
}
  \item{alpha}{numeric. 1 minus the confidence level used in calculating
    the c.i.s 
}
  \item{sim.res}{Character. What simulation samples should be
    returned. If \code{'none'} (the default) the function returns a list
    of 3 arrays (see under 'value'). If \code{'rates'} it returns an
    array of dimension \code{nrow(nd)} x \code{length(mod)} x \code{nB}
    of bootstrap samples of the rates. If \code{'crisk'} it returns an
    array of dimension \code{(nrow(nd)+1)} x \code{length(mod)} x
    \code{nB} of bootstrap samples of the culmulative rates.  Only the
    first letter matters, regardless of whether it is in upper lower
    case.
}  
\item{probs}{Three-way array of simulated cumulative risks classified by
  1) time points, 2) causes (incl. surv) and 3) Samples. A structure as
  returned by \code{ci.Crisk} with \code{sim.res='crisk'}.
}
}
\value{A named list of three-way arrays with results from simulation
  (parametric bootstrap) from the distribution of the parameters in the
  models in \code{mods}:
    
  \itemize{
    \item \code{Crisk} Cumulative risks for the \code{length(mods)}
    events \emph{and} the survival
    \item \code{Srisk} Stacked versions of the cumulative risks
    \item \code{Stime} Sojourn times in each states
  }
  All three arrays have (almost) the same dimensions:
  \itemize{
    \item \code{time}: end points of intervals starting with
    "\code{0}". Length \code{nrow(nd)+1}, except for \code{Stime}
    where it is only \code{nrow(nd)}, "\code{0}" not included.
    
    \item State. \code{Crisk} and \code{Stime} has values
    \code{Surv} plus the names of the list \code{mods} (first
    argument). \code{Srisk} has length \code{length(mod)}, with each
    level representing a cumultive sum of cumulatieve risks, in order
    indicated by the \code{perm} argument. 

    \item \code{ci.50\%}, \code{ci.2.5\%}, \code{ci.97.5\%} representing
    quantiles of the quantities derived from the bootstrap
    samples. If \code{alpha} is different from 0.05, names are of
    course different.
  }
}
\author{
  Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
  \code{\link{mat2pol}}
  \code{\link{simLexis}}
  \code{\link{plotCIF}}
  \code{\link{ci.surv}}
}
\examples{
library(Epi)
data(DMlate)

# A Lexis object for survival
Ldm <- Lexis(entry = list( per = dodm,
                           age = dodm-dobth, 
                           tfd = 0 ),
              exit = list( per = dox ),
       exit.status = factor( !is.na(dodth), labels = c("DM","Dead") ),
              data = DMlate[sample(1:nrow(DMlate),1000),] )
summary(Ldm, timeScales = TRUE)

# Cut at OAD and Ins times
Mdm <- mcutLexis( Ldm,
                   wh = c('dooad','doins'),
           new.states = c('OAD','Ins'),
            precursor = 'Alive',
           seq.states = FALSE,
                 ties = TRUE )
summary( Mdm$lex.dur )

# restrict to DM state
Sdm <- splitLexis(factorize(subset(Mdm, lex.Cst == "DM")),
                  time.scale = "tfd",
                  breaks = seq(0,20,1/12))
summary(Sdm)
summary(Relevel(Sdm, c(1, 4, 2, 3)))

boxes(Relevel(Sdm, c(1, 4, 2, 3)), 
      boxpos = list(x = c(15, 85, 80, 15),
                    y = c(85, 85, 20, 15)),
      scale.R = 100)

# glm models for the cause-specific rates
system.time(
mD <- glm.Lexis(Sdm, ~ Ns(tfd, knots=0:6*2), to = 'Dead') )
system.time(
mO <- glm.Lexis(Sdm, ~ Ns(tfd, knots=0:6*2), to = 'OAD' ) )
system.time(
mI <- glm.Lexis(Sdm, ~ Ns(tfd, knots=0:6*2), to = 'Ins' ) )

# intervals for calculation of predicted rates
int <- 1/100
nd <- data.frame( tfd = seq(int,10,int)-int/2 ) # not the same as the split, 
                                                # and totally unrelated to it

# cumulaive risks with confidence intervals
# (too few timepoints, too few simluations)
system.time(
res <- ci.Crisk(list(OAD = mO, 
                     Ins = mI, 
                    Dead = mD),
                            nd = data.frame(tfd = (1:100-0.5)/10),
                            nB = 100,
                          perm = 4:1))
str(res)
}
\keyword{regression}
