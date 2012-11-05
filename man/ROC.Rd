\name{ROC}
\alias{ROC}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function to compute and draw ROC-curves.}
\description{
Computes sensitivity, specificity and positive and negative predictive
values for a test based on dichotomizing along the variable
\code{test}, for prediction of \code{stat}. Alternatively a
model formula may given, in which case the the linear predictor is the
test variable and the response is taken as the true status variable.
Plots curves of these and a ROC-curve.
}
\usage{
ROC( test = NULL,
     stat = NULL,
     form = NULL,
     plot = c("sp", "ROC"),
       PS = is.null(test),    
       PV = TRUE,             
       MX = TRUE,             
       MI = TRUE,             
      AUC = TRUE,             
     grid = seq(0,100,10),
 col.grid = gray( 0.9 ),
     cuts = NULL,
      lwd = 2,
     data = parent.frame(), 
      ... )
}
\arguments{
  \item{test}{ Numerical variable used for prediction. }
  \item{stat}{ Logical variable of true status. }
  \item{form}{ Formula used in a logistic regression. If this is given,
               \code{test} and \code{stat} are ignored. If not given then
               both \code{test} and \code{stat} must be supplied. }
  \item{plot}{ Character variable. If "sp", the a plot of sensitivity,
specificity and predictive values against test is produced, if "ROC" a
ROC-curve is plotted. Both may be given.}
  \item{PS}{logical, if TRUE the x-axis in the
    plot "ps"-plot is the the predicted probability for
    \code{stat}==TRUE, otherwise it is the scale of \code{test} if this
    is given otherwise the scale of the linear predictor from the
    logistic regression.}
  \item{PV}{Should sensitivity, specificity and
    predictive values at the optimal cutpoint be given on the ROC plot? }
  \item{MX}{Should the ``optimal cutpoint'' (i.e. where sens+spec is
maximal) be indicated on the ROC curve?}
  \item{MI}{Should model summary from the logistic
    regression model be printed in the plot?}
  \item{AUC}{Should the area under the curve (AUC) be printed in the ROC
    plot?}
  \item{grid}{Numeric or logical. If FALSE no background grid is
    drawn. Otherwise a grid is drawn on both axes at \code{grid} percent.}
  \item{col.grid}{Colour of the grid lines drawn.}
  \item{cuts}{Points on the test-scale to be annotated on the
    ROC-curve. }
  \item{lwd}{Thickness of the curves}
  \item{data}{Data frame in which to interpret the variables.}
  \item{\dots}{Additional arguments for the plotting of the
    ROC-curve. Passed on to \code{plot}}
}
%\details{
%}
\value{
  A list with two components:
  \item{res}{dataframe with variables sn, sp, pvp, pvn and fv. The
latter is the unique values of test (for PS==FALSE ) or linear predictor from the
logistic regression}
  \item{lr}{glm object with the logistic regression result used for
construction of the ROC curve}
  0, 1 or 2 plots are produced according to the setting of \code{plot}.
}
\author{Bendix Carstensen, Steno Diabetes Center \& University of
  Copenhagen,
  \url{http://www.biostat.ku.dk/~bxc} }
\examples{
x <- rnorm( 100 )
z <- rnorm( 100 )
w <- rnorm( 100 )
tigol <- function( x ) 1 - ( 1 + exp( x ) )^(-1)
y <- rbinom( 100, 1, tigol( 0.3 + 3*x + 5*z + 7*w ) )
ROC( form = y ~ x + z, plot="ROC" )
}
\keyword{manip}
\keyword{htest}
%\keyword{ROC-curves}
%\keyword{sensitivity}
%\keyword{specificity}
%\keyword{predictive values}

















