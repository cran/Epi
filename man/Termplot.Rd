\name{Termplot}
\alias{Termplot}
\title{
  A wrapper for \code{termplot} that optionally (but by default)
  exponentiates terms, and plot them on a common log-scale. Also scales
  x-axes to the same physical scale.
}
\description{
  The function uses \code{\link{termplot}} to extract terms from a model
  with, say, spline, terms, including the standard errors, computes
  confidence intervals and transform these to the rate / rate-ratio
  scale. Thus the default use is for models on the log-scale such as
  Poisson-regression models. The function produces a plot with panels
  side-by-side, one panel per term, and returns the 
}
\usage{
  Termplot( obj,
           plot = TRUE,
           xlab = NULL,
           ylab = NULL,
            xeq = TRUE,
           yshr = 1,
          alpha = 0.05,
          terms = NULL,
         max.pt = NULL )
}
\arguments{
  \item{obj}{An object with a \code{terms}-method --- for details the
    the documentation for \code{\link{termplot}}. }
  \item{plot}{Should a plot be produced?}
  \item{xlab}{Labels for the \code{x}-axes. Defaults to the names of the terms.}
  \item{ylab}{Labels for the \code{x}-axes. Defaults to blank.}
  \item{xeq}{Should the units all all plots have the same physical scale
    for the \code{x}-axes).}
  \item{yshr}{Shrinking of \code{y}-axis. By default, the \code{y}-axes
    have an extent that accommodates the entire range of confidence
    intervals. This is a shrinking parameter for the \code{y}-axes,
    setting it to less than 1 will lose a bit of the confidence limits
    on some of the panels.}
  \item{alpha}{1 minus the confidence level for computing confidence
    intervals}
  \item{terms}{Which terms should be reported. Passed on to
    \code{\link{termplot}} and eventually \code{\link{predict}}.}
  \item{max.pt}{The maximal number of points in which to report the
    terms. If \code{NULL} all unique points from the analysis dataset
    are reported for each term (this is a feature of \code{\link{termplot}}).}
}
\value{
  A list with one component per term in the model object \code{obj},
  each component is a 4-column matrix with $x$ as the first column, and
  3 columns with estimae and lower and upper confidence limit.
}
\author{
  Bendix Cartensen
}
\seealso{
 \code{\link{Ns}}, \code{termplot}
}
\examples{
# Get the diabetes data and set up as Lexis object
data(DMlate)
DMlate <- DMlate[sample(1:nrow(DMlate),500),]
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit = list(Per=dox),
        exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
               data = DMlate )

# Split in 1-year age intervals
dms <- splitLexis( dml, time.scale="Age", breaks=0:100 )

# Model with 6 knots for both age and period
n.kn <- 6
# Model age-specific rates with period referenced to 2004
( a.kn <- with( subset(dms,lex.Xst=="Dead"),
                quantile( Age+lex.dur, probs=(1:n.kn-0.5)/n.kn ) ) )
( p.kn <- with( subset(dms,lex.Xst=="Dead"),
                quantile( Per+lex.dur, probs=(1:n.kn-0.5)/n.kn ) ) )
m2 <- glm( lex.Xst=="Dead" ~ -1 +
                             Ns( Age, kn=a.kn, intercept=TRUE ) +
                             Ns( Per, kn=p.kn, ref=2004 ),
           offset = log( lex.dur ), family=poisson, data=dms )

# Finally we can plot the two effects:
Termplot( m2, yshr=0.9 )
}
\keyword{hplot}
