\name{Icens}
\alias{Icens}
\alias{print.Icens}
\alias{summary.Icens}
\title{
  Fits a regression model to interval censored data.
  }
\description{   
  The models fitted assumes a piecewise constant baseline rate in
  intervals specified by the argument \code{breaks}, and for the
  covariates either a multiplicative relative risk function (default) or
  an additive excess risk function.
  }
\usage{
  Icens(first.well, last.well, first.ill,
        formula, model.type = c("MRR", "AER"), breaks,
        boot = FALSE, alpha = 0.05, keep.sample = FALSE,
        data)
\method{summary}{Icens}(object, scale = 1, ...)
\method{print}{Icens}(x, scale = 1, digits = 4, ...)
  }
\arguments{
  \item{first.well}{Time of entry to the study, i.e. the time first seen
    without event. Numerical vector.}
  \item{last.well}{Time last seen without event. Numerical vector.}
  \item{first.ill}{Time first seen with event. Numerical vector.}
  \item{formula}{Model formula for the log rate ratios / excess rates.}
  \item{model.type}{Which model should be fitted.}
  \item{breaks}{Breakpoints between intervals in which the underlying
    timescale effect is assumed constant. Any observation outside the range of
    \code{breaks} is discarded.}
  \item{boot}{Should bootstrap be performed to produce confidence
    intervals for parameters. If a number is given this will be the
    number of bootstrap samples. The default is 1000.}
  \item{alpha}{1 minus the confidence level.}
  \item{keep.sample}{Should the bootstrap sample of the parameter values
    be returned?} 
  \item{data}{Data frame in which the times and formula are
    interpreted.}
  \item{object}{an \code{Icens} object.}
  \item{x}{an \code{Icens} object.}
  \item{scale}{scaling factor for rates.}
  \item{digits}{how many digits is used for printing results.}
  \item{ ... }{Other parameters passed on.}  
  }
\details{
  The model is fitted by calling either \code{\link{fit.mult}} or
  \code{\link{fit.add}}. 
  }
\value{
  An object of class \code{"Icens"}: a list with two or three components:
  \item{rates}{A glm object from a binomial model with log-link,
    estimating the baseline rates, and the excess rates if \code{"AER"}
    is specified. The parameters in the linear predictor are rates,
    \emph{not} log-rates. They are in units corresponding to
    (the inverse of) the units
    in which the arguments \code{first.well},  \code{last.well},
    \code{first.ill} and  \code{breaks} are given.}
  \item{cov}{A glm object from a binomial model with complementary
    log-log link, estimating the log-rate-ratios. Only if \code{"MRR"}
    is specified.}
  \item{niter}{Number of iterations, a scalar}
  \item{boot.ci}{If \code{boot=TRUE}, a 3-column matrix with estimates
    and 1-\code{alpha} confidence intervals for the parameters in the model.}
  \item{sample}{A matrix of the parameter estimates from the
    bootstrapping. Rows refer to parameters, columns to bootstrap samples.}
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
  \code{\link{fit.add}}
  \code{\link{fit.mult}}
  }
\examples{
data(hivDK)
# Convert the dates to fractional years so that rates are
# expressed in cases per year and not per day
for(i in 2:4) hivDK[,i] <- cal.yr(hivDK[,i])

# multiplicative rate ratio model, baseline rates and RRs
m.RR <- Icens(entry, well, ill,
              model = "MRR",
            formula = ~ pyr + us,
             breaks = seq(1980, 1990, 5),
               data = hivDK)
# the MRR model returns a list with 2 glm objects.
round(ci.lin(m.RR$rates), 4)         # rates per 1 year
round(ci.lin(m.RR$cov, Exp=TRUE), 4) # RRs
# There is a print method that does the job:
print(m.RR)
# and scales only the rate parameters, here per 10 years
print(m.RR, scale = 10)

# additive excess rate model, baseline rates and excess rates
m.ER <- Icens(entry, well, ill,
              model = "AER",
            formula = ~ pyr + us,
             breaks = seq(1980, 1990, 5),
               data = hivDK)
# The print method returns rates per 1 year:
print(m.ER)
# or per 10 years:
print(m.ER, scale = 10)
  }
\keyword{models}
\keyword{regression}
\keyword{survival}
