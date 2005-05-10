\name{ftrend}
\alias{ftrend}
\title{Fit a floating trend to a factor in generalized linear model}
\description{
  Fits a "floating trend" model to the given factor in a glm in a generalized
  linear model by centering covariates.
}
\usage{
ftrend(object, ...)
}
\arguments{
  \item{object}{fitted \code{lm} or \code{glm} object. The model must not have an intercept term}
  \item{...}{arguments to the \code{nlm} function}
}
\details{
  \code{ftrend()} calculates "floating trend" estimates for factors in
  generalized linear models.  This is an alternative to treatment
  contrasts suggested by Greenland et al. (1999).  If a regression model
  is fitted with no intercept term, then contrasts are  not used for the
  first factor in the model. Instead, there is one parameter  for each
  level of this factor.  However, the interpretation of these
  parameters,  and their variance-covariance matrix, depends on the
  numerical coding used for the  covariates. If an arbitrary constant is
  added to the covariate values, then   the variance matrix is changed. 

  The \code{ftrend()} function takes the fitted model and works out an optimal 
  constant to add to the covariate values so that the covariance matrix is
  approximately diagonal.  The parameter estimates can then be treated as
  approximately independent, thus simplifying their presentation. This is
  particularly useful for graphical display of dose-response relationships
  (hence the name).

  Greenland et al. (1999) originally suggested centring the covariates so that
  their weighted mean, using the fitted weights from the model, is zero.  This
  heuristic criterion is improved upon by \code{ftrend()} which uses the same 
  minimum information divergence criterion as used by Plummer (2003) for
  floating variance calculations. \code{ftrend()} calls \code{nlm()} to
  do the minimization and will pass optional arguments to control it.
}
\note{
  The "floating trend" method is an alternative to the "floating
  absolute risk" method, which is implemented in the function
  \code{float()}.
}
\value{
  A list with the following components
  \item{coef}{coefficients for model with adjusted covariates.}
  \item{vcov}{Variance-covariance matrix of adjusted coefficients.}
}
\references{
  Greenland S, Michels KB, Robins JM, Poole C and Willet WC (1999)
  Presenting statistical uncertainty in trends and dose-response relations,
  \emph{American Journal of Epidemiology}, \bold{149}, 1077-1086.
}
\author{Martyn Plummer}
\seealso{\code{\link{float}}}
\keyword{regression}
