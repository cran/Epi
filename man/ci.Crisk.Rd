\name{ci.Crisk}
\alias{ci.Crisk}
\title{
Compute cumulative risks and expected sojourn times from models for
cause-specific rates.  
}
\description{
  
Consider a list of parametric models for rates of competing events, such
as different causes of death, A, B, C, say. From estimates of the
cause-specific rates we can compute 1) the cumulative risk of being in
each state ('Surv' (=no event) and A, B and C) at different times, 2)
the stacked cumulative rates such as A, A+C, A+C+Surv and 3) the
expected (truncated) sojourn times in each state up to each time point.

This can be done by simple numerical integration using estimates from
models for the cause specific rates. But the standard errors of the
results are analytically intractable.

The function \code{ci.Crisk} computes estimates with confidence
intervals using simulated samples from the parameter vectors of supplied
model objects. Some call this a parametric bootstrap.

The times and other covariates determining the cause-specific rates must
be supplied in a data frame which will be used for predicting rates for
all transitions.
}
\usage{
ci.Crisk(mods,
           nd,
         tnam = names(nd)[1],
           nB = 1000,
         perm = length(mods):0 + 1,
        alpha = 0.05, 
      sim.res = 'none')
}
\arguments{
  \item{mods}{A named list of \code{glm}/\code{gam} model objects
    representing the cause-specific rates. If the list is not named the
    function will crash. The names will be used as names for the states
    (competing risks), while the state without any event will be called
    "\code{Surv}".
  }
  \item{nd}{A data frame of prediction points and covariates to be used
    on all models supplied in \code{mods}. 
  }
  \item{tnam}{Name of the column in \code{nd} which is the time scale.It
    must represent endpoints of equidistant intervals.
  }
  \item{nB}{Scalar. The number of simulations from the (posterior)
    distribution of the model parameters to be used in computing
    confidence limits.  
  }
  \item{perm}{Numerical vector of length \code{length(mods)+1} indicating
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
    array of dimension \code{nrow(nd)} x \code{length(mod)+1} x
    \code{nB} of bootstrap samples of the cumulative rates.  Only the
    first letter matters, regardless of whether it is in upper lower
    case.
}  
}

\value{If \code{sim.res='none'} a named list with 4 components, the first 3
  are 3-way arrays classified by time, state and estimate/confidence
  interval:  
   
  \itemize{
    \item \code{Crisk} Cumulative risks for the \code{length(mods)}
    events \emph{and} the survival
    \item \code{Srisk} Stacked versions of the cumulative risks
    \item \code{Stime} Sojourn times in each states
    \item \code{time} Endpoints of intervals. It is just the numerical
    version of the names of the first dimension of the three arrays
  }
  All three arrays have (almost) the same dimensions:
  \itemize{
    \item time, named as \code{tnam}; endpoints of intervals. Length
    \code{nrow(nd)}.
    
    \item \code{cause}. The arrays \code{Crisk} and \code{Stime} have
    values "\code{Surv}" plus the names of the list \code{mods} (first
    argument). \code{Srisk} has length \code{length(mod)}, with each
    level representing a cumulative sum of cumulative risks, in order
    indicated by the \code{perm} argument.

    \item Unnamed, \code{ci.50\%}, \code{ci.2.5\%}, \code{ci.97.5\%}
    representing quantiles of the quantities derived from the bootstrap
    samples. If \code{alpha} is different from 0.05, names are of course
    different.
  }
  If \code{sim.res='rates'} the function returns bootstrap samples of
    rates for each cause as an array
    classified by time, cause and bootstrap sample.
  
  If \code{sim.res='crisk'} the function returns bootstrap samples of
    cumulative risks for each cause (including survival) as an array
    classified by time, state (= causes + surv) and bootstrap sample.
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
Mdm <- mcutLexis(Ldm,
                  wh = c('dooad','doins'),
          new.states = c('OAD','Ins'),
          seq.states = FALSE,
                ties = TRUE)
summary(Mdm$lex.dur)

# restrict to DM state and split
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
int <- 1 / 100
nd <- data.frame(tfd = seq(0, 10, int)) # not the same as the split, 
                                        # and totally unrelated to it

# cumulaive risks with confidence intervals
# (too few timepoints, too few simluations)
system.time(
res <- ci.Crisk(list(OAD = mO, 
                     Ins = mI, 
                    Dead = mD),
                            nd = data.frame(tfd = 0:100 / 10),
                            nB = 100,
                          perm = 4:1))
str(res)
}
\keyword{regression}
