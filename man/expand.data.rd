\name{expand.data}
\alias{expand.data}
\title{
  Function to expand data for regression analysis of interval censored
  data.
  }
\description{
  This is a utility function.
  
  The original records with \code{first.well}, \code{last.well} and
  \code{first.ill} are
  expanded to multiple records; several for each interval where the
  person is known to be well and one where the person is known to fail.
  At the same time columns for the covariates needed to estimate rates
  and the response variable are generated.
  }
\usage{
expand.data(fu, formula, breaks, data)
}
\arguments{
  \item{fu}{A 3-column matrix with \code{first.well}, \code{last.well} and
  \code{first.ill} in each row.}
  \item{formula}{Model fromula, used to derive the model matrix.}
  \item{breaks}{Defines the intervals in which the baseline rate is
    assumed constant. All follow-up before the first and after the
    last break is discarded.}
  \item{data}{Datafrem in which \code{fu} and \code{formula} is interpreted.}
}
\value{
  Returns a list with three components
  \item{rates.frame}{Dataframe of covariates for estimation of the
    baseline rates --- one per interval defined by \code{breaks}.}
  \item{cov.frame}{Dataframe for estimation of the covariate effects. A
    data-framed version of the designmatrix from \code{formula}.}
  \item{y}{Response vector.}
}
\references{
  B Carstensen: Regression models for interval censored
  survival data: application to HIV infection in Danish homosexual
  men. Statistics in Medicine, 15(20):2177-2189, 1996.
  }
\author{
  Martyn Plummer, \email{plummer@iarc.fr}
  }
\seealso{
  \code{\link{Icens}}
  \code{\link{fit.mult}}
  \code{\link{fit.add}}
  }
\examples{
  }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
