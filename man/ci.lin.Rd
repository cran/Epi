\name{ci.lin}
\alias{ci.lin}
\title{
  Compute linear functions of parameters with s.e.
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
      diffs = FALSE,
       fnam = !diffs,
       vcov = FALSE,
      alpha = 0.05,
        Exp = FALSE )
}
\arguments{
  \item{obj}{A model object (of class \code{lm}, \code{glm}, \code{lme},
    \code{coxph} or \code{polr}). 
  }
  \item{ctr.mat}{Contrast matrix to be multiplied to the parameter
    vector, i.e. the desired linear function of the parameters.}
  \item{subset}{The subset of the parameters to be used. If given as a
    character vector, the elements are in turn matched against the
    parameter names (using \code{grep}) to find the subset. Repeat
    parameters may result from using a character vector. This is
    considered a facility.} 
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
    returned? If this is set, \code{Exp} is ignored.}
  \item{alpha}{Significance level for the confidence intervals.}
  \item{Exp}{If \code{TRUE} columns 5:6 are replaced with exp( columns 1,5,6 ).}
}
\value{
  A matrix with number of rows and rownames as \code{ctr.mat}. The
  columns are Estimate, Std.Err, z, P, 2.5\% and 97.5\%.
  If \code{vcov=TRUE} a list with components \code{est}, the desired
  functional of the parameters and \code{vcov}, the variance
  covariance matrix of this, is returned but not printed.
  If code{Exp==TRUE} the confidence intervals for the parameters are
  replaced with three columns: exp(estimate,c.i.).
}
\author{
  Bendix Carstensen,
  \url{http://www.pubhealth.ku.dk/~bxc}
}
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
}
\keyword{models}
\keyword{regression}
