\name{ci.eta}
\alias{ci.eta}
\title{
  Linear predictor (\code{eta}) from a formula, coefficients, vcov and a
  prediction frame.
}
\description{
  Computes the linear predictor with its confidence limits from the
  model formula and the estimated parameters with the vcov.
}
\usage{
ci.eta(form, cf, vcv, newdata,
       alpha = 0.05, df = Inf, raw = FALSE)
}
\arguments{
  \item{form}{A model formula. A one-sided formula will suffice; left
    side will be ignored if two-sided.
}
  \item{cf}{Coefficients from a model using \code{formula}.
}
  \item{vcv}{variance-covariance matrix from a model using \code{formula}.
}
  \item{newdata}{Prediction data frame with variables used in
    \code{formula}. Can also be a list of 2 or 4 prediction frames, for
    details see \code{\link{ci.lin}}.
}
  \item{alpha}{Significance level for calculation of c.i.
}
  \item{df}{Integer. Number of degrees of freedom in the t-distribution
          used to compute the quantiles used to construct the
          confidence intervals.
}
  \item{raw}{Logical. Should predictions and their vcov be returned
    instead of predictions and confidence limits?
}
}
\details{
  Does pretty much the same as \code{\link{ci.lin}}, but requires only a
  formula and coefficients with vcov. Designed to avoid saving entire
  (homongously large) model objects and still be able to compute
  predictions. But only the linear predictor is returned, if there is a
  link in your model function it is your own responsibility to
  back-transform. If the model formula contains reference to vectors of
  spline knots or similar these must be in the global environment.

  There is no guarantee that this function works for models that do not
  inherit from \code{lm}. But there is a guantee that it will not work
  for \code{gam} objects.
}
\value{
The linear predictor for the \code{newdata} with a confidence interval
as a 3 column matrix. If \code{raw=TRUE} the linear predictor and its
variance-covariance matrix.
}
\author{
  Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{\code{\link{ci.lin}}
}
\keyword{models}
\keyword{regression}
\concept{prediction}
\concept{prediction frame}
\concept{linear predictor}
