\name{fit.add}
\alias{fit.add}
\title{
  Fit an additive excess risk model to interval censored data.
  }
\description{
  Utility function.
    
  The model fitted assumes a piecewise constant intensity for the
  baseline, and that the covariates act additively on the rate scale.
  }
\usage{
  fit.add( y, rates.frame, cov.frame, start )
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
\value{
  A list with one component:
  \item{rates}{A glm object from a binomial model with log-link function.}
}
\references{
  B Carstensen: Regression models for interval censored
  survival data: application to HIV infection in Danish homosexual
  men. Statistics in Medicine, 15(20):2177-2189, 1996.

  CP Farrington: Interval censored survival data: a generalized linear
  modelling approach. Statistics in Medicine, 15(3):283-292, 1996.
  }
\author{
  Martyn Plummer, \email{martyn.plummer@r-project.org}
  }
\seealso{
  \code{\link{Icens}}
  \code{\link{fit.mult}}
  }
\examples{
  data( HIV.dk ) 
  }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
