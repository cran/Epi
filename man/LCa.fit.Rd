\name{LCa.fit}
\alias{LCa.fit}
\alias{print.LCa}
\alias{summary.LCa}
\alias{plot.LCa}
\alias{predict.LCa}
\title{
  Fit Lee-Carter-type models for rates to arbitrarily shaped observations
  of rates in a Lexis diagram. 
}
\description{
  The Lee-Carter model is originally defined as a model for rates
  observed in A-sets (age by period) of a Lexis diagram, as
  log(rate(x,t)) = a(x) + b(x)k(t), using one parameter per age(x) and
  period(t). This function uses natural splines for a(), b() and k(),
  placing knots for each effect such that the number of events is the
  same between knots. 
}
\usage{
LCa.fit( data, A, P, D, Y,
         model = "APa",    # or one of "ACa", "APaC", "APCa" or "APaCa" 
         a.ref,            # age reference for the interactions
        pi.ref = a.ref,    # age reference for the period interaction
        ci.ref = a.ref,    # age reference for the cohort interaction
         p.ref,            # period reference for the interaction
         c.ref,            # cohort reference for the interactions
          npar = c(a = 6,  # no. knots for main age-effect
                   p = 6,  # no. knots for period-effect
                   c = 6,  # no. knots for cohort-effect
                  pi = 6,  # no. knots for age in the period interaction
                  ci = 6), # no. knots for age in the cohort interaction
            VC = TRUE,     # numerical calculation of the Hessian?
         alpha = 0.05,     # 1 minus confidence level
           eps = 1e-6,     # convergence criterion
         maxit = 100,      # max. no iterations
         quiet = TRUE )    # cut the crap
\method{print}{LCa}( x, ... )
\method{summary}{LCa}( object, show.est=FALSE, ... )
\method{plot}{LCa}( x, ... )
\method{predict}{LCa}( object, newdata,
                        alpha = 0.05,
                        level = 1-alpha,
                          sim = ( "vcov" \%in\% names(object) ),
                          ... )
      }
\arguments{
  \item{data}{
  A data frame. Must have columns \code{A}(age), \code{P}(period, that is
  calendar time), \code{D}(no. of events) and \code{Y}(person-time,
  exposure). Alternatively these four quantities can be given as
  separate vectors:
  }
  \item{A}{
  Vector of ages (midpoint of observation).
  }
  \item{P}{
  Vector of period (midpoint of observation).
  }
  \item{D}{
  Vector of no. of events.
  }
  \item{Y}{
  Vector of person-time. Demographers would say "exposure", bewildering epidemiologists.
  }
  \item{a.ref}{
  Reference age for the age-interaction term(s) \code{pi(x)} and/or
  \code{pi(x)}, where \code{pi(a.ref)=1} and \code{ci(a.ref)=1}.
  }
  \item{pi.ref}{
  Same, but specifically for the interaction with period.
  }
  \item{ci.ref}{
  Same, but specifically for the interaction with cohort.
  }
  \item{p.ref}{
  Reference period for the time-interaction term \code{kp(t)} where \code{kp(p.ref)=0}.
  }
  \item{c.ref}{
  Reference period for the time-interaction term \code{kp(t)} where \code{kc(c.ref)=0}.
  }
  \item{model}{
  Character, either \code{"APa"} which is the classical Lee-Carter model
  for log-rates, other possibilities are \code{"ACa"}, \code{"APCa"},
  \code{"APaC"} or \code{"APaCa"}, see details.
  }
  \item{npar}{
  A (possibly named) vector or list, with either the number of knots or
  the actual vectors of knots for each term. If unnamed, components are
  taken to be in the order (a,b,t), if the model is "APaCa" in the order
  (a,p,c,pi,ci).  If a vector, the three integers indicate the number of
  knots for each term; these will be placed so that there is an equal
  number of events (\code{D}) between each, and half as many below the
  first and above the last knot. If \code{npar} is a list of scalars the
  behavior is the same. If \code{npar} is a list of vectors, these are
  taken as the knots for the natural splines. See details for naming
  convention. 
  }
  \item{VC}{
  Logical. Should the variance-covariance matrix of the parameters be
  computed by numerical differentiation? See details.
  }
  \item{alpha}{
  1 minus the confidence level used when calculating
  confidence intervals for estimates in \code{LCa.fit} and for
  predictions by \code{predict.LCa}.
  }
  \item{eps}{
  Convergence criterion for the deviance, we use the the relative
  difference between deviance from the two models fitted.
  }
  \item{maxit}{
  Maximal number of iterations.
  }
  \item{quiet}{
  Shall I shut up or talk extensively to you about iteration progression etc.?
  }
  \item{object}{An \code{LCa} object, see under "Value".}
  \item{show.est}{Logical. Should the estimates be printed?}
  \item{x}{An \code{LCa} object, see under "Value".}
  \item{newdata}{Prediction data frame, must have columns \code{A} and
                 \code{P}. Any \code{Y} column is ignored, predictions are given in
  units of the \code{Y} supplied for the call that generated the
  \code{LCa} object.}
  \item{level}{Confidence level.}
  \item{sim}{Logical or numeric. If \code{TRUE}, prediction c.i.s will be
  based on 1000 simulations from the posterior parameters. If numeric,
  it will be based on that number of simulations.}
  \item{...}{Additional parameters passed on to the method.}
}
\details{
  The Lee-Carter model is non-linear in age and time so does not fit
  in the classical glm-Poisson framework. But for fixed \code{b(x)} it
  is a glm, and also for fixed \code{a(x)}, \code{k(t)}. The function
  alternately fits the two versions until the same fit is produced (same
  deviance).

  The multiplicative age by period term could equally well have been a
  multiplicative age by cohort or even both. Thus the most extensive
  model has 5 continuous functions:

  \deqn{\log(\lambda(a,p)) = f(a) + b_p(a)k_p(p) + b_c(a)k_c(p-a)}{%
         log( lambda(a,p)) = f(a) + b_p(a)k_p(p) + b_c(a)k_c(p-a)}

  Each of these is fitted by a natural spline, with knots placed at the
  quantiles of the events along the age (a), calendar time (p) respective
  cohort (p-a) scales. Alternatively the knots can be specified explicitly
  in the argument \code{npar} as a named list, where
  \code{a} refers to \eqn{f(a)}{f(a)},
  \code{p} refers to \eqn{k_p(p)}{k_p(p)},
  \code{c} refers to \eqn{k_c(p-a)}{k_c(p-a)},
  \code{pi} (\code{p}eriod \code{i}nteraction) refers to \eqn{b_p(a)}{b_p(a)}
  and
  \code{ci} (\code{c}ohort \code{i}nteraction) refers to \eqn{b_c(p-a)}{b_c(p-a)}.

  The naming convention for the models is a capital \code{P} and/or
  \code{C} if the effect is in the model followed by a lower case
  \code{a} if there is an interaction with age. Thus there are 5 different
  models that can be fitted: \code{APa}, \code{ACa}, \code{APaC}  \code{APCa}
  and \code{APaCa}.
	  
  The standard errors of the parameters from the two separate model fits
  in the iterations are however wrong; they are conditional on a subset
  of the parameters having a fixed value. However, analytic calculation
  of the Hessian is a bit of a nightmare, so this is done numerically
  using the \code{hessian} function from the \code{numDeriv} package if
  \code{VC=TRUE}.
  
  The coefficients and the variance-covariance matrix of these are used
  in \code{predict.LCa} for a parametric bootstrap (that is, a
  simulation from a multivariate normal with mean equal to the parameter
  estimates and variance as the estimated variance-covariance) to get
  confidence intervals for the predictions if \code{sim} is \code{TRUE}
  --- which it is by default if they are part of the object.

  The \code{plot} for \code{LCa} objects merely produces between 3 and 5
  panels showing each of the terms in the model. These are mainly for
  preliminary inspection; real reporting of the effects should use
  proper relative scaling of the effects.}

\value{
  \code{LCa.fit} returns an object of class \code{LCa} (smooth
  effects \code{L}ee-\code{Ca}rter model); it is a list with the
  following components:
\item{model}{Character, either \code{APa}, \code{ACa}, \code{APaC},
  \code{APCa} or \code{APaCa}, indicating the variable(s) interacting
  with age.} 
\item{ax}{3-column matrix of age-effects, c.i. from the age-time
  model. Row names are the unique occurring ages in the
  dataset. Estimates are rates.}  
\item{pi}{3-column matrix of age-period interaction effects, c.i. from the age
  model. Row names are the actually occurring ages in the
  dataset. Estimates are multipliers of the log-RRs in \code{kp},
  centered at 1 at \code{pi.ref}.}
\item{kp}{3-column matrix of period-effects, with c.i.s from the
  age-time model. Row names are the actually occurring times in the 
  dataset. Estimates are rate-ratios centered at 1 at \code{p.ref}.}
\item{ci}{3-column matrix of age-cohort interaction effects, c.i. from the age
  model. Row names are the actually occurring ages in the
  dataset. Estimates are multipliers of the log-RRs in \code{kc},
  centered at 1 at \code{ci.ref}.}
\item{kc}{3-column matrix of cohort-effects, with c.i.s from the age-time
  model. Row names are the actually occurring times in the
  dataset. Estimates are rate-ratios centered at 1 at \code{c.ref}.}
\item{mod.at}{\code{glm} object with the final age-time model --- estimates
  the terms \code{ax}, \code{kp}, \code{kc}. Gives
  the same fit as the \code{mod.b} model after convergence.}
\item{mod.b}{\code{glm} object with the final age model --- estimates
  the terms \code{pi}, \code{ci}. Gives
  the same fit as the \code{mod.at} model after convergence.}
\item{coef}{All coefficients from both models, in the order \code{ax},
  \code{kp}, \code{kc}, \code{pi}, \code{ci}. Only present if
  \code{LCa.fit} were called with \code{VC=TRUE} (the default).} 
\item{vcov}{Variance-covariance matrix of coefficients from both
  models, in the same order as in the \code{coef}. Only present if
  \code{LCa.fit} were called with \code{VC=TRUE}.} 
\item{knots}{List of vectors of knots used in for the age, period and
  cohort effects.}
\item{refs}{List of reference points used for the age, period and
  cohort terms in the interactions.}
\item{deviance}{Deviance of the model}
\item{df.residual}{Residual degrees of freedom}
\item{iter}{Number of iterations used to reach convergence.}

  \code{plot.LCa} plots the estimated effects in separate panels,
  using a log-scale for the baseline rates (\code{ax}) and the time-RR
  (\code{kt}). For the \code{APaCa} model 5 panels are plotted.
  
  \code{summary.LCa} returns (invisibly) a matrix with the parameters
  from the models and a column of the conditional se.s and additionally
  of the se.s derived from the numerically computed Hessian (if
  \code{LCa.fit} were called with \code{VC=TRUE}.)

  \code{predict.LCa} returns a matrix with one row per row in
  \code{newdata}. If \code{LCa.fit} were called with \code{VC=TRUE}
  there will be 3 columns, namely prediction (1st column) and c.i.s
  based on a simulation of parameters from a multivariate normal with
  mean \code{coef} and variance \code{vcov} using the median and
  \code{alpha}/2 quantiles from the \code{sim} simulations.  If
  \code{LCa.fit} were called with \code{VC=FALSE} there will be 6
  columns, namely estimates and c.i.s from age-time model
  (\code{mod.at}), and from the age-interaction model (\code{mod.b}),
  both using conditional variances, and therefore likely with too narrow
  confidence limits.
  }
\author{
  Bendix Carstensen, \url{http://bendixcarstensen.com}

  This function was conceived while teaching a course on APC models at
  the Max Planck Institute of Demographic Research (MPIDR,
  \url{https://www.demogr.mpg.de/en/}) in Rostock in May 2016
  (\url{http://bendixcarstensen.com/APC/MPIDR-2016/}), and finished
  during a week long research stay there, kindly sponsored by the MPIDR.
  }
\seealso{
 \code{\link{apc.fit}},
 \code{\link{apc.LCa}},
 \code{\link[demography]{lca}} }
\examples{
library( Epi )
# Load the testis cancer data by Lexis triangles
data( testisDK )
tc <- subset( testisDK, A>14 & A<60 )
head( tc )

# We want to see rates per 100,000 PY
tc$Y <- tc$Y / 10^5

# Fit the Lee-Carter model with age-period interaction (default)
LCa.tc <- LCa.fit( tc, model="ACa", a.ref=30, p.ref=1980, quiet=FALSE, eps=10e-4, maxit=50 )

LCa.tc
summary( LCa.tc )

# Inspect what we got
names( LCa.tc )

# show the estimated effects
par( mfrow=c(1,3) )
plot( LCa.tc )

# Prediction data frame for ages 15 to 60 for two time points: 
nd <- data.frame( A=15:60 )
# LCa predictions
p70 <- predict.LCa( LCa.tc, newdata=cbind(nd,P=1970), sim=1000 )
p90 <- predict.LCa( LCa.tc, newdata=cbind(nd,P=1990), sim=1000 )

# Inspect the curves from the parametric bootstrap (simulation):
par( mfrow=c(1,1) )
head( cbind(p70,p90) )
matplot( nd$A, cbind(p70,p90), type="l", lwd=c(6,3,3), lty=c(1,3,3),
         col=rep( 2:3, each=3 ), log="y",
         ylab="Testis cancer incidence per 100,000 PY in 1970 resp. 1990", xlab="Age" )
}
\keyword{models}
\keyword{regression}
