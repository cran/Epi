\name{clogistic}
\alias{clogistic}
\title{Conditional logistic regression}
\description{
  Estimates a logistic regression model by maximizing the conditional
  likelihood.  The conditional likelihood calculations are exact, and
  scale efficiently to strata with large numbers of cases.
}
\usage{
clogistic(formula, strata, data, subset, na.action, init,
model = TRUE, x = FALSE, y = TRUE, contrasts = NULL,
iter.max=20, eps=1e-6, toler.chol = sqrt(.Machine$double.eps)) 
}
\arguments{
  \item{formula}{Model formula}
  \item{strata}{Factor describing membership of strata for conditioning}
  \item{data}{data frame containing the variables in the formula and
    strata arguments}
  \item{subset}{subset of records to use}
  \item{na.action}{missing value handling}
  \item{init}{initial values}
  \item{model}{
    a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value}
  \item{x,y}{
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.
  }
  \item{contrasts}{
    an optional list. See the \code{contrasts.arg} of
    \code{model.matrix.default}
  }
  \item{iter.max}{maximum number of iterations}
  \item{eps}{
    Convergence tolerance.  Iteration continues until the relative
    change in the conditional log likelihood is less than \code{eps}.
    Must be positive.
  }
  \item{toler.chol}{
    Tolerance used for detection of a singularity during a Cholesky
    decomposition of the variance matrix. This is used to detect
    redundant predictor variables. Must be less than \code{eps}.
  }
}
\value{
  An object of class \code{"clogistic"}. This is a list containing
  the following components:
  \item{coefficients}{
    the estimates of the log-odds ratio parameters.  If the model is
    over-determined there will be missing values in the vector corresponding
    to the redundant columns in the model matrix.
  }
  \item{var}{
    the variance matrix of the coefficients.  Rows and columns corresponding to 
    any missing coefficients are set to zero. 
  }
  \item{loglik}{
    a vector of length 2 containing the log-likelihood with the initial
    values and with the final values of the coefficients. 
  }
  \item{iter}{
    number of iterations used. 
  }
  \item{n}{
    number of observations used. Observations may be dropped
    either because they are missing, or because they belong to a
    homogeneous stratum. For more details on which observations were
    used, see \code{informative} below.
  }
  \item{informative}{
    if \code{model=TRUE}, a logical vector of length equal to the number
    of rows in the model frame.  This indicates whether an observation
    is informative. Strata that are homogeneous with respect to either
    the outcome variable or the predictor variables are uninformative,
    in the sense that they can be removed without modifying the
    estimates or standard errors. If \code{model=FALSE}, this is NULL.
  }
  
  The output will also contain the following, for documentation see the
  \code{glm} object: \code{terms}, \code{formula},
  \code{call},  \code{contrasts}, \code{xlevels}, and, optionally,
  \code{x}, \code{y}, and/or \code{frame}. 
}
\examples{
  data(bdendo)
  clogistic(d ~ cest + dur, strata=set, data=bdendo)
}
\author{Martyn Plummer}
\seealso{\code{\link{glm}}}
\keyword{models}
