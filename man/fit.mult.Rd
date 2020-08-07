\name{fit.mult}
\alias{fit.mult}
\title{
  Fits a multiplicative relative risk model to interval censored data.
  }
\description{
  Utility function.
  
  The model fitted assumes a piecewise constant baseline rate in
  intervals specified by the argument \code{breaks}, and a
  multiplicative relative risk function.  
}
\usage{
  fit.mult( y, rates.frame, cov.frame, start )
  }
\arguments{
  \item{y}{Binary vector of outcomes}
  \item{rates.frame}{Dataframe expanded from the original data by
    \code{\link{expand.data}}, cooresponding to covariates for the rate
    parameters.}
  \item{cov.frame}{ do., but covariates corresponding to the
    \code{formula} argument of \code{\link{Icens}}} 
  \item{start}{Starting values for the rate parameters. If not supplied,
       then starting values are generated.}
  }
\details{
  The model is fitted by alternating between two generalized linear
  models where one estimates the underlying rates in the intervals, and
  the other estimates the log-relative risks.
  }
\value{
  A list with three components:
  \item{rates}{A glm object from a binomial model with log-link,
    estimating the baseline rates.}
  \item{cov}{A glm object from a binomial model with complementary
    log-log link, estimating the log-rate-ratios}
  \item{niter}{Nuber of iterations, a scalar}
  }
\references{
  B Carstensen: Regression models for interval censored
  survival data: application to HIV infection in Danish homosexual
  men. Statistics in Medicine, 15(20):2177-2189, 1996.

  CP Farrington: Interval censored survival data: a generalized linear
  modelling approach. Statistics in Medicine, 15(3):283-292, 1996.
  }
\author{
  Martyn Plummer, \email{martyn.plummer@r-project.org},
  Bendix Carstensen, \email{b@bxc.dk} }
\seealso{
  \code{\link{Icens}}
  \code{\link{fit.add}}
  }
\examples{
  data( HIV.dk ) 
  }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
