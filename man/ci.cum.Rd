\name{ci.cum}
\alias{ci.cum}
\title{ Compute cumulative sum of estimates. }
\description{
Computes the cumulative sum of parameter functions and the
standard error of it. Optionally the exponential is applied to the
parameter functions before it is cumulated.
}
\usage{
ci.cum( obj,
    ctr.mat = NULL,
     subset = NULL,
       intl = 1,
      alpha = 0.05,
        Exp = TRUE,
     sample = FALSE )
}
\arguments{
  \item{obj}{A model object (of class
             \code{lm}, \code{glm},
             \code{coxph}, \code{survreg},
             \code{lme},\code{mer},\code{nls},\code{gnlm},
             \code{MIresult}
          or \code{polr}). }
  \item{ctr.mat}{ Contrast matrix defining the parameter functions from
    the parameters of the model. }
  \item{subset}{ Subset of the parameters of the model to which
    \code{ctr.mat} should be applied. }
  \item{intl}{ Interval length for the cumulation. Either a constant or
    a numerical vector of length \code{nrow(ctr.mat)}. }
  \item{alpha}{ Significance level used when computing confidence limits. }
  \item{Exp}{ Should the parameter function be exponentiated before it is
    cumulated?}
  \item{sample}{Should a sample of the original parameters be used to
    compute a cumulative rate?}
}
\details{
The purpose of this function is to compute cumulative rate based on a
model for the rates. If the model is a multiplicative model for the
rates, the purpose of \code{ctr.mat} is to return a vector of rates or
log-rates when applied to the coefficients of the model. If log-rates
are returned from the model, the they should be exponentiated before
cumulated, and the variances computed accordingly. Since log-linear
models are the most common the \code{Exp} parameter defaults to TRUE.
}
\value{
   A matrix with 4 columns: Estimate, lower and upper c.i. and standard
   error. If \code{sample} is TRUE, a sampled vector is reurned, if
   \code{sample} is numeric a matrix with \code{sample} columns is
   returned, each column a cumulative rate based on a random sample from
   the distribution of the parameter estimates.
 }
\author{
  Bendix Carstensen,
  \url{http://www.pubhealth.ku.dk/~bxc}
}
\seealso{ See also \code{\link{ci.lin}} }
\examples{
# Packages required for this example
library( splines )
library( survival )
data( lung )
par( mfrow=c(1,2) )

# Plot the Kaplan-meier-estimator
plot( survfit( Surv( time, status==2 ) ~ 1, data=lung ) )

# Declare data as Lexis
lungL <- Lexis( exit=list("tfd"=time),
                exit.status=(status==2)*1, data=lung )
summary( lungL )

# Cut the follow-up every 10 days
sL <- splitLexis( lungL, "tfd", breaks=seq(0,1100,10) )
str( sL )
summary( sL )

# Fit a Poisson model with a natural spline for the effect of time.
# Extract the variables needed
D <- status(sL, "exit")
Y <- dur(sL)
tB <- timeBand( sL, "tfd", "left" )
MM <- ns( tB, knots=c(50,100,200,400,700), intercept=TRUE )
mp <- glm( D ~ MM - 1 + offset(log(Y)),
           family=poisson, eps=10^-8, maxit=25 )

# Contrast matrix to extract effects, i.e. matrix to multiply with the
# coefficients to produce the log-rates: unique rows of MM, in time order.
T.pt <- sort( unique( tB ) )
T.wh <- match( T.pt, tB )
Lambda <- ci.cum( mp, ctr.mat=MM[T.wh,], intl=diff(c(0,T.pt)) )

# Put the estimated survival function on top of the KM-estimator
matlines( c(0,T.pt[-1]), exp(-Lambda[,1:3]), lwd=c(3,1,1), lty=1, col="Red" )

# Extract and plot the fitted intensity function
lambda <- ci.lin( mp, ctr.mat=MM[T.wh,], Exp=TRUE )
matplot( T.pt, lambda[,5:7]*10^3, type="l", lwd=c(3,1,1), col="black", lty=1,
         log="y", ylim=c(0.2,20) )
}
\keyword{models}
\keyword{regression}
