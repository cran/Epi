\name{plotEst}
\alias{plotEst}
\alias{pointsEst}
\alias{linesEst}
\title{
  Plot estimates with confidence limits (forest plot)
}
\description{
  Plots parameter estimates with confidence intervals, annotated with
  parameter names. A dot is plotted at the estimate and a horizontal
  line extending from the lower to the upper limit is superimposed.
}
\usage{
plotEst( ests,
            y = dim(ests)[1]:1,
          txt = rownames(ests),
       txtpos = y,
         ylim = range(y)-c(0.5,0),
         xlab = "",
         xtic = nice(ests[!is.na(ests)], log = xlog),
         xlim = range( xtic ),
         xlog = FALSE,
          pch = 16,
          cex = 1,
          lwd = 2,
          col = "black",
      col.txt = "black",
     font.txt = 1,
    col.lines = col,
   col.points = col,
         vref = NULL,
         grid = FALSE,
     col.grid = gray(0.9),
  restore.par = TRUE,
          ... )

linesEst( ests, y = dim(ests)[1]:1, pch = 16, cex = 1, lwd = 2,
          col="black", col.lines=col, col.points=col, ... )

pointsEst( ests, y = dim(ests)[1]:1, pch = 16, cex = 1, lwd = 2,
          col="black", col.lines=col, col.points=col, ... )
}
\arguments{
  \item{ests}{Matrix with three columns: Estimate, lower limit, upper
    limit. If a model object is supplied, \code{\link{ci.lin}} is
    invoked for this object first.}
  \item{y}{Vertical position of the lines.}
  \item{txt}{Annotation of the estimates. Either a character vector or
          an expression vector.}
  \item{txtpos}{Vertical position of the text. Defaults to \code{y}.}
  \item{ylim}{Extent of the vertical axis.}
  \item{xlab}{Annotation of the horizontal axis.}
  \item{xtic}{Location of tickmarks on the x-axis.}
  \item{xlim}{Extent of the x-axis.}
  \item{xlog}{Should the x-axis be logarithmic?}
  \item{pch}{What symbol should be used?}
  \item{cex}{Expansion of the symbol.}
  \item{col}{Colour of the points and lines.}
  \item{col.txt}{Colour of the text annotating the estimates.}
  \item{font.txt}{Font for the text annotating the estimates.}
  \item{col.lines}{Colour of the lines.}
  \item{col.points}{Colour of the symbol.}
  \item{lwd}{Thickness of the lines.}
  \item{vref}{Where should vertical reference line(s) be drawn?}
  \item{grid}{If TRUE, vertical gridlines are drawn at the
    tickmarks. If a numerical vector is given vertical lines are drawn
    at \code{grid}.}
  \item{col.grid}{Colour of the vertical gridlines}
  \item{restore.par}{Should the graphics parameters be restored? If set
    to \code{FALSE} the coordinate system will still be available for
    additional plotting, and \code{par("mai")} will still have the very
    large value set in order to make room for the labelling of the
    estimates.}
  \item{...}{Arguments passed on to \code{ci.lin} when a model object is
             supplied as \code{ests}.}
}
\details{
  \code{plotEst} makes a news plot, whereas \code{linesEst} and
  \code{pointsEst} (identical functions) adds to an existing plot.

  If a model object of class \code{"glm"}, \code{"coxph"}, \code{"clogistic"} or
  \code{"gnlm"} is supplied the argument \code{xlog} defaults to
  \code{TRUE}, and exponentiated estimates are extracted by default.
}
\value{
  NULL
}
\author{
  Bendix Carstensen,
  \email{b@bxc.dk},
  \url{http://bendixcarstensen.com}}
\seealso{
  ci.lin
}
\examples{
# Bogus data and a linear model
f <- factor( sample( letters[1:5], 100, replace=TRUE ) )
x <- rnorm( 100 )
y <- 5 + 2 * as.integer( f ) + 0.8 * x + rnorm(100) * 2
m1 <- lm( y ~ f )

# Produce some confidence intervals for contrast to first level
( cf <- ci.lin( m1, subset=-1 )[,-(2:4)] )

# Plots with increasing amounts of bells and whistles
par( mfcol=c(3,2), mar=c(3,3,2,1) )
plotEst( cf )
plotEst( cf, grid=TRUE, cex=2, lwd=3 )
plotEst( cf, grid=TRUE, cex=2, col.points="red", col.lines="green" )
plotEst( cf, grid=TRUE, cex=2, col.points="red", col.lines="green",
             xlog=TRUE, xtic=c(1:8), xlim=c(0.8,6) )
rownames( cf )[1] <- "Contrast to fa:\n fb"
plotEst( cf, grid=TRUE, cex=2, col.points=rainbow(4),
                                col.lines=rainbow(4), vref=1 )
#
etxt <- expression("Plain text, qouted",
                   "combined with maths:"*sqrt(a)*phi[c],
                   f^d*"  Hb"*A[1][c],
                   eff^e*"  kg/"*m^2)
plotEst( cf, txt=etxt, grid=TRUE, cex=2, col.points=rainbow(4),
                                         col.lines =rainbow(4), vref=1 )
  }
\keyword{hplot}
\keyword{models}
