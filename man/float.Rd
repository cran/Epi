\name{float}
\alias{float}
\alias{print.floated}
\title{Calculate floated variances}
\description{
  Given a fitted  model object, the \code{float()} function calculates
  floating variances (aka quasi-variances) for a given factor in the model.
}
\usage{
float(object, factor, iter.max=50)
}
\arguments{
  \item{object}{a fitted model object}
  \item{factor}{character string giving the name of the factor of
    interest. If this is not given, the first factor in the model is used.}
  \item{iter.max}{Maximum number of iterations for EM algorithm} }
\details{
  float() implements the "floating absolute risk" proposal of Easton,
  Peto and Babiker(1992). This is an alternative way of presenting
  parameter estimates for factors in regression models, which avoids
  some of the difficulties of treatment contrasts. It was originally
  designed for epidemiological studies of relative risk (hence the name)
  but the idea is widely applicable.  

  A problem with treatment contrasts is that they are not
  orthogonal. The variances of the treatment contrasts may be inflated by
  a poor choice of reference level, and the correlations between them
  may be very high.  float() associates each level of the factor,
  including the reference level, with a"floating" variance (or
  quasi-variance). Floating variances are not real variances, but
  they can be used to calculate the variance of any contrast by treating
  each level as independent.
  
  Plummer (2003) showed that floating variances can be derived from a
  covariance structure model applied to the variance-covariance
  matrix of the parameter estimates. This model can be fitted by
  minimizing the Kullback-Leibler information divergence between the
  true and distributions for the parameter estimates and the
  distribution given by the covariance structure model. Fitting is
  done using the EM algorithm.

  In order to check the goodness-of-fit of the floating variance
  model, \code{float()} compares the standard errors predicted
  by the model with the standard errors derived from the true
  variance-covariance matrix of the parameter contrasts. The maximum
  and minimum ratios between true and model standard errors are
  calculated over all possible contrasts. These should be within 5
  percent, or the use of the floating variances may lead to invalid
  confidence intervals.
}
\value{
  An object of class \code{floated}. This is a list with the following
  components
  \item{coef}{vector of coefficients. These are the same as the
    treatment contrasts but the reference level is present with
    coefficient 0.}
  \item{var}{vector of floating (or quasi-) variances}
  \item{limits}{Bounds on the accuracy of standard errors over all
    possible contrasts}
}
\note{
  Menezes(1999) and Firth and Menezes (2004) take a slightly different
  approach to this problem, using a pseudo-likelihood approach to fit
  the quasi-variance model. Their work is implemented in the package qvcalc. 
}
\references{
  Easton DF, Peto J and Babiker GAG (1991) Floating absolute risk: An
  alternative to relative risk in survival and case control analysis
  avoiding an arbitrary reference group. \emph{Statistics in Medicine},
  \bold{10}, 1025-1035.

  Firth D and Mezezes RX (2004)  Quasi-variances.
  \emph{Biometrika} \bold{91}, 65-80.

  Menezes RX(1999)  More useful standard errors for group and factor
  effects in generalized linear models.  \emph{D.Phil. Thesis},
  Department of Statistics, University of Oxford.

  Plummer M (2003) Improved estimates of floating absolute risk,
  \emph{Statistics in Medicine}, \bold{23}, 93-104.
}
\author{Martyn Plummer}
\seealso{\code{\link{ftrend}}, \code{qvcalc}}
\keyword{regression}
