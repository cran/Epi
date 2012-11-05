\name{fit.baseline}
\alias{fit.baseline}
\title{
  Fit a piecewise contsnt intesity model for interval censored data.
  }
\description{
  Utility function
  
  Fits a binomial model with logaritmic link, with \code{y} as outcome
  and covariates in \code{rates.frame} to estimate rates in the
  inttervals between \code{breaks}.
  }
\usage{
fit.baseline( y, rates.frame, start )
}
\arguments{
  \item{y}{Binary vector of outcomes}
  \item{rates.frame}{Dataframe expanded from the original data by
    \code{\link{expand.data}}} 
  \item{start}{Starting values for the rate parameters. If not supplied,
       then starting values are generated.}
}
\value{
  A \code{\link{glm}} object, with binomial error and logaritmic link.
  }
\author{
  Martyn Plummer, \email{plummer@iarc.fr}
  }
\seealso{
  \code{\link{fit.add}}
  \code{\link{fit.mult}}
  }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
