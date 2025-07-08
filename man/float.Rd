\name{float}
\alias{float}
\alias{print.floated}
\title{Calculate floated variances}
\description{
  Given a fitted  model object, the \code{float()} function calculates
  floating variances (a.k.a. quasi-variances) for a given factor in the model.
}
\usage{
float(object, factor, iter.max=50, ...)
}
\arguments{
  \item{object}{a fitted model object.}
  \item{factor}{character string giving the name of the factor of
    interest. If this is not given, the first factor in the model is used.}
  \item{iter.max}{Maximum number of iterations for EM algorithm.}
  \item{...}{Optional arguments passed to the \code{vcov()} method for
  the fitted model object.}
}
\details{

  The \code{float()} function implements the "floating absolute risk"
  proposal of Easton, Peto and Babiker (1992). This is an alternative way
  of presenting parameter estimates for factors in regression models,
  which avoids some of the difficulties of treatment contrasts. It was
  originally designed for epidemiological studies of relative risk, but
  the idea is widely applicable.

  Treatment contrasts are not orthogonal. Consequently, the variances of
  treatment contrast estimates may be inflated by a poor choice of
  reference level, and the correlations between them may also be high.
  The \code{float()} function associates each level of the factor with a
  "floating" variance (or quasi-variance), including the reference
  level.  Floating variances are not real variances, but they can be
  used to calculate the variance error of contrast by treating each
  level as independent.

  Plummer (2003) showed that floating variances can be derived from a
  covariance structure model applied to the variance-covariance matrix
  of the contrast estimates. This model can be fitted by minimizing the
  Kullback-Leibler information divergence between the true distribution
  of the parameter estimates and the simplified distribution given by
  the covariance structure model. Fitting is done using the EM
  algorithm.

  In order to check the goodness-of-fit of the floating variance model,
  the \code{float()} function compares the standard errors predicted by
  the model with the standard errors derived from the true
  variance-covariance matrix of the parameter contrasts. The maximum and
  minimum ratios between true and model-based standard errors are
  calculated over all possible contrasts. These should be within 5
  percent, or the use of the floating variances may lead to invalid
  confidence intervals.

}
\value{
  An object of class \code{floated}. This is a list with the following
  components
  \item{coef}{A vector of coefficients. These are the same as the
    treatment contrasts but the reference level is present with
    coefficient 0.}
  \item{var}{A vector of floating (or quasi-) variances}
  \item{limits}{The bounds on the accuracy of standard errors over all
    possible contrasts}
}
\note{
  Menezes(1999) and Firth and Menezes (2004) take a slightly different
  approach to this problem, using a pseudo-likelihood approach to fit
  the quasi-variance model. Their work is implemented in the package
  qvcalc.
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
