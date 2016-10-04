\name{ci.lin}
\alias{ci.lin}
\alias{ci.mat}
\alias{ci.exp}
\alias{ci.pred}
\alias{Wald}
\title{
  Compute linear functions of parameters with standard errors and
  confidence limits
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
      diffs = FALSE,
       fnam = !diffs,
       vcov = FALSE,
      alpha = 0.05,
         df = Inf,
        Exp = FALSE,
     sample = FALSE )
ci.exp( ..., Exp = TRUE, pval=FALSE )
Wald( obj, H0=0, ... )
ci.mat( alpha = 0.05, df = Inf )
ci.pred( obj, newdata,
         Exp = NULL,
       alpha = 0.05,
          df = Inf )
}
\arguments{
  \item{obj}{A model object (of class
    \code{lm},
    \code{glm},
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
  \item{ctr.mat}{Contrast matrix to be multiplied to the parameter
    vector, i.e. the desired linear function of the parameters.}
  \item{subset}{The subset of the parameters to be used. If given as a
    character vector, the elements are in turn matched against the
    parameter names (using \code{grep}) to find the subset. Repeat
    parameters may result from using a character vector. This is
    considered a facility.}
  \item{subint}{\code{sub}set selection like for \code{subset}, except that
    elements of a character vector given as argument will be used to
    select a number of subsets of parameters and only the \code{int}ersection
    of these is returned.}
  \item{diffs}{If TRUE, all differences between parameters
    in the subset are computed. \code{ctr.mat} is ignored. If \code{obj}
    inherits from \code{lm}, and \code{subset} is given as a string
    \code{subset} is used to search among the factors in the model and
    differences of all factor levels for the first match are shown.
    If \code{subset} does not match any of the factors in the model, all
    pairwise differences between parameters matching are returned.}
  \item{fnam}{Should the common part of the parameter names be included
    with the annotation of contrasts? Ignored if \code{diffs==T}. If a
    sting is supplied this will be prefixed to the labels.}
  \item{vcov}{Should the covariance matrix of the set of parameters be
    returned? If this is set, \code{Exp} is ignored. See details.}
  \item{alpha}{Significance level for the confidence intervals.}
  \item{df}{Integer. Number of degrees of freedom in the t-distribution used
    to compute the quantiles used to construct the confidence intervals.}
  \item{Exp}{If \code{TRUE} columns 5:6 are replaced with exp( columns
    1,5,6 ). For \code{ci.pred} it indicates whether the predictions
    should be exponentiated - the default is to make a prediction on the
    scale of the linear predictor and transform it by the inverse link
    function; if FALSE, the prediction on the link scale is returned.}
  \item{sample}{Logical or numerical. If \code{TRUE} or numerical a
    sample of size \code{as.numeric(sample)} is drawn from the
    multivariate normal with mean equal to the (\code{subet} defined)
    coefficients and variance equal to the estimated variance-covriance
    of these. These are then transformed by \code{ctr.mat} and
    returned.}
  \item{pval}{Logical. Should a column of P-values be included with the estimates
    and confidence intervals output by \code{ci.exp}.}
  \item{H0}{Numeric. The null values for the selected/transformed parameters to
    be tested by a Wald test. Must have the same length as the selected
    parameter vector.}
  \item{\ldots}{Parameters passed on to \code{ci.lin}.}
  \item{newdata}{Data frame of covariates where prediction is made.}
  }
\value{
  \code{ci.lin} returns a matrix with number of rows and row names as
  \code{ctr.mat}. The columns are Estimate, Std.Err, z, P, 2.5\% and
  97.5\%.  If \code{vcov=TRUE} a list with components \code{est}, the
  desired functional of the parameters and \code{vcov}, the variance
  covariance matrix of this, is returned but not printed.  If
  \code{Exp==TRUE} the confidence intervals for the parameters are
  replaced with three columns: exp(estimate,c.i.).

  \code{ci.exp} returns only the exponentiated parameter estimates with
  confidence intervals. It is merely a wrapper for \code{ci.lin},
  fishing out the last 3 columns from \code{ci.lin(...,Exp=TRUE)}. If
  you just want the estimates and confidence limits, but not
  exponentiated, use \code{ci.exp(...,Exp=FALSE)}.

  \code{Wald} computes a Wald test for a subset of (possibly linear
  combinations of) parameters being equal to the vector of null
  values. The selection of the subset of parameters is the same as for
  \code{ci.lin}. Using the \code{ctr.mat} argument makes it possible to
  do a Wald test for equality of parameters. \code{Wald} returns a named
  numerical vector of length 3, with names \code{Chisq}, \code{d.f.} and
  \code{P}.

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
  }
\author{
  Bendix Carstensen,
  \url{BendixCarstensen.com} &
  Michael Hills
  }
\seealso{See also \code{\link{ci.cum}} for a function computing
  cumulative sums of (functions of) parameter estimates.}
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

# A Wald test of whether the g-parameters are 0
Wald( mm, subset="g" )
# Wald test of whether the three first f-parameters are equal:
( CM <- rbind( c(1,-1,0,0), c(1,0,-1,0)) )
Wald( mm, subset="f", ctr.mat=CM )
# or alternatively
( CM <- rbind( c(1,-1,0,0), c(0,1,-1,0)) )
Wald( mm, subset="f", ctr.mat=CM )
}
\keyword{models}
\keyword{regression}
