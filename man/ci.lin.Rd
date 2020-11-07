\name{ci.lin}
\alias{ci.lin}
\alias{ci.mat}
\alias{ci.exp}
\alias{ci.pred}
\alias{ci.ratio}
\alias{Wald}
\title{
  Compute linear functions of parameters with standard errors and
  confidence limits, optionally transforming to a different scale.
  }
\description{
  For a given model object the function computes a linear function of
  the parameters and the corresponding standard errors, p-values and
  confidence intervals.
  }
\usage{
ci.lin( obj,
    ctr.mat = NULL,
     subset = NULL,
     subint = NULL,
      xvars = NULL,
      diffs = FALSE,
       fnam = !diffs,
       vcov = FALSE,
      alpha = 0.05,
         df = Inf,
        Exp = FALSE,
     sample = FALSE )
ci.exp( ..., Exp = TRUE, pval = FALSE )
Wald( obj, H0=0, ... )
ci.mat( alpha = 0.05, df = Inf )
ci.pred( obj, newdata,
         Exp = NULL,
       alpha = 0.05 )
ci.ratio( r1, r2,
         se1 = NULL,
         se2 = NULL,
      log.tr = !is.null(se1) & !is.null(se2),
       alpha = 0.05,
        pval = FALSE )
}
\arguments{
  \item{obj}{A model object (in general of class \code{glm}, but for
    \code{ci.lin} and \code{ci.exp} it may also be of class
    \code{lm},
    \code{coxph},
    \code{survreg},
    \code{clogistic},
    \code{cch},
    \code{lme},
    \code{mer},
    \code{lmerMod},
    \code{gls},
    \code{nls},
    \code{gnlm},
    \code{MIresult},
    \code{mipo},
    \code{polr},
 or \code{rq}).
  }
  \item{ctr.mat}{Matrix, data frame or list (of two or four data frames).

    If \code{ctr.mat} is a matrix, it should be a contrast matrix to be
    multiplied to the parameter vector, i.e. the desired linear function
    of the parameters.

    If it is a data frame it should have columns corresponding to a
    prediction frame, see details.
    
    If it is a list, it must contain two or four data frames that are
    (possibly partial) prediction frames for \code{obj}, see argument
    \code{xvars} and details.}
  \item{xvars}{Character vector. If quantitative variables in the model are omitted
    from data frames supplied in a list to \code{ctr.mat}, they should
    be listed here. Omitted factors need not be mentioned here.}
  \item{subset}{The subset of the parameters to be used. If given as a
    character vector, the elements are in turn matched against the
    parameter names (using \code{grep}) to find the subset. Repeat
    parameters may result from using a character vector. This is
    considered a facility.}
  \item{subint}{Character. \code{sub}set selection, but where each
    element of the character vector is used to
    select a subset of parameters and only the \code{int}ersection
    of these is returned.}
  \item{diffs}{If TRUE, all differences between parameters
    in the subset are computed, and the \code{subset} argumnt is
    required. The argument \code{ctr.mat} is ignored. If \code{obj}
    inherits from \code{lm}, and \code{subset} is given as a string
    \code{subset} is used to search among the factors in the model and
    differences of all factor levels for the first match are shown.
    If \code{subset} does not match any of the factors in the model, all
    pairwise differences between parameters matching are returned.}
  \item{fnam}{Should the common part of the parameter names be included
    with the annotation of contrasts? Ignored if \code{diffs==T}. If a
    string is supplied this will be prefixed to the labels.}
  \item{vcov}{Should the covariance matrix of the set of parameters be
    returned? If this is set, \code{Exp} is ignored. See details.}
  \item{alpha}{Significance level for the confidence intervals.}
  \item{df}{Integer. Number of degrees of freedom in the t-distribution used
    to compute the quantiles used to construct the confidence intervals.}
  \item{Exp}{For \code{ci.lin}, if \code{TRUE} columns 5:6 are replaced
    with exp( columns 1,5,6 ). For \code{ci.exp}, if \code{FALSE}, the
    untransformed parameters are returned. For \code{ci.pred} it
    indicates whether the predictions should be exponentiated - the
    default (\code{Exp=NULL}) is to make a prediction with a Wald CI on
    the scale of the linear predictor and back-transform it by the
    inverse link function; if \code{FALSE}, the prediction on the link
    scale is returned.}
  \item{sample}{Logical or numerical. If \code{TRUE} or numerical a
    sample of size \code{as.numeric(sample)} is drawn from the
    multivariate normal with mean equal to the (\code{subset} defined)
    coefficients and variance equal to the estimated variance-covariance
    of these. These are then transformed by \code{ctr.mat} and
    returned.}
  \item{pval}{Logical. Should a column of P-values be included with the
    estimates and confidence intervals output by \code{ci.exp}.}
  \item{H0}{Numeric. The null values for the selected/transformed
    parameters to be tested by a Wald test. Must have the same length as
    the selected parameter vector.}
  \item{\ldots}{Parameters passed on to \code{ci.lin}.}
  \item{newdata}{Data frame of covariates where prediction is made.}
  \item{r1,r2}{Estimates of rates in two independent groups, with
    confidence limits. Can be either 3-column matrices or data frames
    with estimates and confindece intervals or 2 two column structures
    with confidence limits. Only the confidence limits}
  \item{se1,se2}{Standard errors of log-rates in the two groups. If
    given, it is assumed that \code{r1} and \code{r2} represent
    log-rates.}
  \item{log.tr}{Logical, if true, it is assumed that \code{r1} and
    \code{r2} represent log-rates with confidence intervals.}
  }
\value{
    
  \code{ci.lin} returns a matrix with number of rows and row names as
  \code{ctr.mat}. The columns are Estimate, Std.Err, z, P, 2.5\% and
  97.5\% (or according to the value of \code{alpha}).  If
  \code{vcov=TRUE} a list of length 2 with components \code{coef} (a
  vector), the desired functional of the parameters and \code{vcov} (a
  square matrix), the variance covariance matrix of this, is returned
  but not printed. If \code{Exp==TRUE} the confidence intervals for the
  parameters are replaced with three columns: exp(estimate,c.i.).

  \code{ci.exp} returns only the exponentiated parameter estimates with
  confidence intervals. It is merely a wrapper for \code{ci.lin},
  fishing out the last 3 columns from \code{ci.lin(...,Exp=TRUE)}. If
  you just want the estimates and confidence limits, but not
  exponentiated, use \code{ci.exp(...,Exp=FALSE)}.
  
  If \code{ctr.mat} is a data frame, the model matrix corresponding to
  this is constructed and supplied. This is only supported for objects
  of class \code{lm}, \code{glm}, \code{gam} and \code{coxph}.

  So the default behaviour will be to produce the same as
  \code{ci.pred}, apparently superfluous. The purpose of this is to
  allow the use of the arguments \code{vcov} that produces the
  variance-covariance matrix of the predictions, and \code{sample} that
  produces a sample of predictions using sampling from the multivariate
  normal with mean equal to parameters and variance equal to the
  hessian.

  If \code{ctr.mat} is a list of two data frames, the difference of the
  predictions from using the first versus the last as newdata arguments
  to predict is computed. Columns that would be identical in the two
  data frames can be omitted (see below), but names of numerical
  variables omitted must be supplied in a character vector
  \code{xvars}. Factors omitted need not be named.

  If the second data frame has only one row, this is replicated to match
  the number of rows in the first. This facility is primarily aimed at
  teasing out RRs that are non-linear functions of a quantitative
  variable without setting up contrast matrices using the same code as
  in the model. Note however if splines are used with computed knots
  stored in a vector such as \code{Ns(x,knots=kk)} then the \code{kk}
  must be available in the (global) environment; it will not be found
  inside the model object. In practical terms it means that if you save
  model objects for later prediction you should save the knots used in
  the spline setups too.

  If \code{ctr.mat} is a list of four data frames, the difference of the
  difference of predictions from using the first and second versus
  difference of predictions from using the third and fourth is computed.
  Simply \code{(pr1-pr2) - (pr3-pr4)} with obvious notation. Useful to
  derive esoteric interaction effects.
  
  Finally, only arguments \code{Exp}, \code{vcov}, \code{alpha} and
  \code{sample} from \code{ci.lin} are honored when \code{ctr.mat} is a
  data frame or a list of two data frames.

  You can leave out variables (columns) from the two data frames that
  would be identical, basically variables not relevant for the
  calculation of the contrast. In many cases \code{ci.lin} (really
  \code{Epi:::ci.dfr}) can figure out the names of the omitted columns,
  but occasionally you will have to supply the names of the omitted
  variables in the \code{xvars} argument. Factors omitted need not be
  listed in \code{xvars}, although no harm is done doing so.

  \code{Wald} computes a Wald test for a subset of (possibly linear
  combinations of) parameters being equal to the vector of null
  values as given by \code{H0}. The selection of the subset of
  parameters is the same as for \code{ci.lin}. Using the \code{ctr.mat}
  argument makes it possible to do a Wald test for equality of
  parameters. \code{Wald} returns a named numerical vector of length 3,
  with names \code{Chisq}, \code{d.f.} and \code{P}.

  \code{ci.mat} returns a 2 by 3 matrix with rows \code{c(1,0,0)} and
  \code{c(0,-1,1)*1.96}, devised to post-multiply to a p by 2 matrix with
  columns of estimates and standard errors, so as to produce a p by 3 matrix
  of estimates and confidence limits. Used internally in \code{ci.lin} and
  \code{ci.cum}.
  The 1.96 is replaced by the appropriate quantile from the normal or
  t-distribution when arguments \code{alpha} and/or \code{df} are given.

  \code{ci.pred} returns a 3-column matrix with estimates and upper and
  lower confidence intervals as columns. This is just a convenience
  wrapper for \code{predict.glm(obj,se.fit=TRUE)} which returns a rather
  unhandy structure. The prediction with c.i. is made in the \code{link}
  scale, and by default transformed by the inverse link, since the most
  common use for this is for multiplicative Poisson or binomial models
  with either log or logit link.

  \code{ci.ratio} returns the rate-ratio of two independent set of
  rates given with confidence intervals or s.e.s. If \code{se1} and
  \code{se2} are given and \code{log.tr=FALSE} it is assumed that
  \code{r1} and \code{r2} are rates and \code{se1} and \code{se2} are
  standard errors of the log-rates.
  }
\author{
  Bendix Carstensen,
  \url{http://bendixcarstensen.com} & Michael Hills
}

\seealso{See also \code{\link{ci.cum}} for a function computing
  cumulative sums of (functions of) parameter estimates, and
  \code{\link{ci.surv}} for a function computing confidence intervals
  for survival functions based on smoothed rates. The example code for
  \code{\link{matshade}} has an application of predicting a rate-ratio
  using a list of two prediction frame in the \code{ctr.mat} argument.}

\examples{
# Bogus data:
f <- factor( sample( letters[1:5], 200, replace=TRUE ) )
g <- factor( sample( letters[1:3], 200, replace=TRUE ) )
x <- rnorm( 200 )
y <- 7 + as.integer( f ) * 3 + 2 * x + 1.7 * rnorm( 200 )

# Fit a simple model:
mm <- lm( y ~ x + f + g )
ci.lin( mm )
ci.lin( mm, subset=3:6, diff=TRUE, fnam=FALSE )
ci.lin( mm, subset=3:6, diff=TRUE, fnam=TRUE )
ci.lin( mm, subset="f", diff=TRUE, fnam="f levels:" )
print( ci.lin( mm, subset="g", diff=TRUE, fnam="gee!:", vcov=TRUE ) )

# Use character defined subset to get ALL contrasts:
ci.lin( mm, subset="f", diff=TRUE )

# Suppose the x-effect differs across levels of g:
mi <- update( mm, . ~ . + g:x )
ci.lin( mi )
# RR a vs. b by x:
nda <- data.frame( x=-3:3, g="a", f="b" )
ndb <- data.frame( x=-3:3, g="b", f="b" )
# 
ci.lin( mi, list(nda,ndb) )
# Same result if f column is omitted because "f" columns are identical
ci.lin( mi, list(nda[,-3],ndb[,-3]) )
# However, crashes if knots in spline is supplied, and non-factor omitted
xk <- -1:1
xi <- c(-0.5,0.5)
ww <- rnorm(200)
mi <- update( mm, . ~ . -x + ww + Ns(x,knots=xk) + g:Ns(x,knots=xi) )
# Will crash 
try( cbind( nda$x, ci.lin( mi, list(nda,ndb) ) ) )
# Must specify num vars (not factors) omitted from nda, ndb
cbind( nda$x, ci.lin( mi, list(nda,ndb), xvars="ww" ) )

# A Wald test of whether the g-parameters are 0
Wald( mm, subset="g" )
# Wald test of whether the three first f-parameters are equal:
( CM <- rbind( c(1,-1,0,0), c(1,0,-1,0)) )
Wald( mm, subset="f", ctr.mat=CM )
# or alternatively
( CM <- rbind( c(1,-1,0,0), c(0,1,-1,0)) )
Wald( mm, subset="f", ctr.mat=CM )

# Confidence intervals for ratio of rates
# Rates and ci supplied, but only the range (lower and upper ci) is used
ci.ratio( cbind(10,8,12.5), cbind(5,4,6.25) )
ci.ratio( cbind(8,12.5), cbind(4,6.25) )

# Beware of the offset when making predictions with ci.pred and ci.exp
\dontrun{
library( mgcv )
data( mortDK )
m.arg  <- glm( dt ~ age , offset=log(risk) , family=poisson, data=mortDK )
m.form <- glm( dt ~ age + offset(log(risk)), family=poisson, data=mortDK )
a.arg  <- gam( dt ~ age , offset=log(risk) , family=poisson, data=mortDK )
a.form <- gam( dt ~ age + offset(log(risk)), family=poisson, data=mortDK )

nd <- data.frame( age=60:65, risk=100 )
round( ci.pred( m.arg , nd ), 4 )
round( ci.pred( m.form, nd ), 4 )
round( ci.exp ( m.arg , nd ), 4 )
round( ci.exp ( m.form, nd ), 4 )
round( ci.pred( a.arg , nd ), 4 )
round( ci.pred( a.form, nd ), 4 )
round( ci.exp ( a.arg , nd ), 4 )
round( ci.exp ( a.form, nd ), 4 )

nd <- data.frame( age=60:65 )
try( ci.pred( m.arg , nd ) )
try( ci.pred( m.form, nd ) )
try( ci.exp ( m.arg , nd ) )
try( ci.exp ( m.form, nd ) )
try( ci.pred( a.arg , nd ) )
try( ci.pred( a.form, nd ) )
try( ci.exp ( a.arg , nd ) )
try( ci.exp ( a.form, nd ) )
}
# The offset may be given as an argument (offset=log(risk))
# or as a term (+offset(log)), and depending on whether we are using a
# glm or a gam Poisson model and whether we use ci.pred or ci.exp to
# predict rates the offset is either used or ignored and either
# required or not; the state of affairs can be summarized as:
#
#                     offset
#                     -------------------------------------
#                     usage                 required?
#                     ------------------    ---------------                      
# function  model     argument   term       argument   term
# ---------------------------------------------------------
# ci.pred   glm       used       used       yes        yes
#           gam       ignored    used       no         yes
#  		      
# ci.exp    glm       ignored    ignored    no         yes
#           gam       ignored    ignored    no         yes
# ---------------------------------------------------------
}
\keyword{models}
\keyword{regression}
