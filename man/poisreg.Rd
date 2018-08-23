\name{poisreg}
\alias{poisreg}

\title{Family Object for Poisson Regression}
\usage{
poisreg(link = "log")
}
\arguments{
  \item{link}{a specification for the model link function.  The
    \code{poisreg} family accepts the links \code{identity},
    \code{log} and \code{inverse}.
  }
}
\description{

  The \code{poisreg} family allows Poisson regression models to be
  fitted using the \code{glm} function.
  
  In a Poisson regression model, we assume that the data arise from a
  Poisson process. We observe D disease events in follow up time Y and
  wish to estimate the incidence rate, which is assumed to be constant
  during the follow-up period for any individual. The incidence rate
  varies between individuals according to the predictor variables and
  the link function in the model specification.
                                                                     
  When using the \code{poisreg} family in the \code{glm} function, the
  response should be specified as a two-column matrix with the first
  column giving the number of events (D) and the second column giving
  the observation time (Y). This is similar to the \code{binomial}
  family for which a two-column outcome can be used representing the
  number of successes and the number of failures.

}
\note{

  When using the log link, Poisson regression can also be carried out
  using the \code{poisson} family by including the log follow-up time
  \code{log(Y)} as an offset. However this approach does not generalize
  to other link functions. The \code{poisreg} family allows more general
  link functions including additive risk models with \code{poisreg(link
    = "identity")}. 

}
\value{
  An object of class \code{"family"}.  See \code{\link[stats]{family}}
  for details.

  The family name, represented by the element \code{"family"} in the
  returned object, is \code{"poisson"} and not \code{"poisreg"}. This is
  necessary to prevent the \code{summary.glm} function from estimating
  an overdispersion parameter (which should be fixed at 1) and therefore
  giving incorrect standard errors for the estimates.

}
\examples{
  ## Estimate incidence rate of diabetes in Denmark (1996-2015) by
  ## age and sex
  data(DMepi)
  DMepi$agegrp <- cut(DMepi$A, seq(from=0, to=100, by=5))
  inc.diab <- glm(cbind(X, Y.nD) ~ 0 + agegrp + sex, family=poisreg,
                  data=DMepi)
  ## The coefficients for agegrp are log incidence rates for men in each
  ## age group. The coefficient for sex is the log of the female:male
  ## incidence rate ratio.
  summary(inc.diab)

  ## Smooth function with non-constant M/F RR:
  requireNamespace("mgcv")
  library( mgcv )
  gam.diab <- gam( cbind(X, Y.nD) ~ s(A,by=sex) + sex,
                   family=poisreg,
                   data=DMepi)
  ## There is no need/use for Y.nD in prediction data frames:
  nM <- data.frame( A=20:90, sex="M" )
  nF <- data.frame( A=20:90, sex="F" )
  ## Rates are returned in units of (1 year)^-1, so we must scale the
  ## rates by hand: 
  matshade( nM$A, cbind( ci.pred(gam.diab,     nM    )*1000,
                         ci.pred(gam.diab,        nF )*1000,
                         ci.exp( gam.diab,list(nM,nF)) ),
            plot=TRUE, col=c("blue","red","black"),
            log="y", xlab="Age", ylab="DM incidence rates per 1000     /     M vs. F RR" )
  abline(h=1)
}
\seealso{
  \code{\link[stats]{glm}}, \code{\link[stats]{family}}.
}
\keyword{models}
