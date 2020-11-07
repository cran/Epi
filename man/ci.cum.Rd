\name{ci.cum}
\alias{ci.cum}
\alias{ci.surv}
\title{ Compute cumulative sum of estimates. }
\description{
Computes the cumulative sum of parameter functions and the
standard error of it. Used for computing the cumulative rate, or the
survival function based on a \code{glm} with parametric baseline.
}
\usage{
ci.cum( obj,
    ctr.mat = NULL,
     subset = NULL,
       intl = 1,
      alpha = 0.05,
        Exp = TRUE,
     ci.Exp = FALSE,
     sample = FALSE )
ci.surv( obj,
    ctr.mat = NULL,
     subset = NULL,
       intl = 1,
      alpha = 0.05,
        Exp = TRUE,
     sample = FALSE )
}
\arguments{
  \item{obj}{A model object (of class
             \code{lm}, \code{glm}. }
  \item{ctr.mat}{Matrix or data frame.

 If \code{ctr.mat} is a matrix, it should be a contrast matrix to be
 multiplied to the parameter vector, i.e. the desired linear
 function of the parameters.

 If it is a data frame it should have columns corresponding to a
 prediction data frame for \code{obj}, see details for \code{\link{ci.lin}}.}
  \item{subset}{ Subset of the parameters of the model to which a matrix
    \code{ctr.mat} should be applied. }
  \item{intl}{ Interval length for the cumulation. Either a constant or
    a numerical vector of length \code{nrow(ctr.mat)}. }
  \item{alpha}{ Significance level used when computing confidence limits. }
  \item{Exp}{ Should the parameter function be exponentiated before it is
    cumulated?}
  \item{ci.Exp}{ Should the confidence limits for the cumulative rate be
  computed on the log-scale, thus ensuring that exp(-cum.rate) is always
  in [0,1]?}
  \item{sample}{Should a sample of the original parameters be used to
    compute a cumulative rate?}
}
\details{  
  The purpose of this function is to the compute cumulative rate
  (integrated intensity) at a set of points based on a model for the
  rates. \code{ctr.mat} is a matrix which, when premultiplied to the
  parameters of the model returns the (log)rates at a set of increasing
  time points. If log-rates are returned from the model, the they should
  be exponentiated before cumulated, and the variances computed
  accordingly. Since the primary use is for log-linear Poisson models
  the \code{Exp} parameter defaults to TRUE.

  Each row in the object supplied via \code{ctr.mat} is assumed to
  represent a midpoint in an interval. \code{ci.cum} will then return
  the cumulative rates at the \emph{end} of these
  intervals. \code{ci.surv} will return the survival probability at the
  \emph{start} of each of these intervals, assuming the the first
  interval starts at 0 - the first row of the result is \code{c(1,1,1)}.
  
  The \code{ci.Exp} argument ensures that the confidence intervals for
  the cumulative rates are always positive, so that exp(-cum.rate) is
  always in [0,1].
}
\value{
   A matrix with 3 columns: Estimate, lower and upper c.i. and standard
   error, unless \code{se=TRUE}, in which case the standard error is
   returned too. If \code{sample} is TRUE, a sampled vector is returned, if
   \code{sample} is numeric a matrix with \code{sample} columns is
   returned, each column a cumulative rate based on a random sample from
   the distribution of the parameter estimates.

   \code{ci.surv} returns a 3 column matrix with estimate, lower and
   upper confidence interval.
}
\author{
  Bendix Carstensen,
  \url{http://bendixcarstensen.com}
}
\seealso{ See also \code{\link{ci.lin}}, \code{\link{ci.pred}} }
\examples{
# Packages required for this example
library( splines )
library( survival )
data( lung )
par( mfrow=c(1,2) )

# Plot the Kaplan-meier-estimator
plot( survfit( Surv( time, status==2 ) ~ 1, data=lung ) )

# Declare data as Lexis
lungL <- Lexis( exit=list(tfd=time),
                exit.status=(status==2)*1,
                data=lung )
summary( lungL )

# Split the follow-up every 10 days
sL <- splitLexis( lungL, "tfd", breaks=seq(0,1100,10) )
summary( sL )

# Fit a Poisson model with a natural spline for the effect of time (left
# end points of intervals are used as covariate)
mp <- glm( cbind(lex.Xst==1,lex.dur) ~ Ns(tfd,knots=c(0,50,100,200,400,700)),
           family=poisreg, data=sL )

# mp is now a model for the rates along the time scale tfd
# prediction data frame for select time points on this time scale
nd <- data.frame( tfd = seq(5,995,10) ) # *midpoints* of intervals
Lambda <- ci.cum ( mp, nd, intl=10 )
surv   <- ci.surv( mp, nd, intl=10 )

# Put the estimated survival function on top of the KM-estimator
# recall the ci.surv return the survival at *start* of intervals
matshade( nd$tfd-5, surv, col="Red", alpha=0.15 )

# Extract and plot the fitted intensity function
lambda <- ci.pred( mp, nd )*365.25 # mortality 
matshade( nd$tfd, lambda, log="y", ylim=c(0.2,5), plot=TRUE,
          xlab="Time since diagnosis", ylab="Mortality per year" )

# same thing works with gam from mgcv
library( mgcv )
mg <- gam( cbind(lex.Xst==1,lex.dur) ~ s(tfd), family=poisreg, data=sL )
matshade( nd$tfd-5, ci.surv( mg, nd, intl=10 ), plot=TRUE,
          xlab="Days since diagnosis", ylab="P(survival)" )
matshade( nd$tfd  , ci.pred( mg, nd )*365.25  , plot=TRUE, log="y",
          xlab="Days since diagnosis", ylab="Mortality per 1 py" )
}
\keyword{models}
\keyword{regression}
