\name{tabplot}
\alias{tabplot}
\title{
  Graphical display of a 2-way contingency table
}
\description{
Entries in a table are plotted as rectangles proportional to the entry
in the table. Width of rectangles are proportional to column totals,
height proportional to entries within each column, hence areas are
proportional to entries in the table. 
}
\usage{
  tabplot( M,
         col,
      border = "black",
         lwd = 2,
     collabs = TRUE,
     rowlabs = NULL,
       equal = FALSE,
         las = 1,
        main = NULL,
    cex.main = 1.0,
       vaxis = FALSE )
}
\arguments{
  \item{M}{Two-way table}
  \item{col}{colors to use for coloring within each column. Defaults to
    a grayscale. May also be a function that takes an integer argument,
    as e.g. \code{rainbow()}.}
  \item{border}{color of borders around rectangles.}
  \item{lwd}{width of the lines around rectangles.}
  \item{collabs}{should colums be labelled.}
  \item{rowlabs}{character "r" or "l": rows labelled on left or right side}
  \item{equal}{should colums be plotted of equal width? If yes a plot
    similar to that obtainable from barplot is the result.}
  \item{las}{how should labelling be rotated.}
  \item{main}{heading for the plot}
  \item{cex.main}{character expansion for the heading.}
  \item{vaxis}{should a vertical axis be drawn. If character it gives
    the side where it is drawn.}
}
\details{
  The function offers a few more facilities for two-way tables than
  \code{\link{mosaicplot}}, but is restricted to two-way tables as
  input. 
}
\value{
  NULL. The function is used for its sideeffects.
}
\author{
  Bendix Carstensen, Steno Diabetes Center \& Dept. of Biostatistics,
  University of Copenhagen
  \email{bxc@steno.dk},
  \url{http://www.pubhealth.ku.dk/~bxc}
}
\seealso{
  See Also \code{\link{barplot}}, \code{\link{plot.table}},
  \code{\link{mosaicplot}}
}
\examples{
b <- sample( letters[1:4], 300, replace=TRUE, prob=c(3,1,2,4)/10 )
a <- rnorm( 300 ) - as.integer( factor( b ) ) / 8
tb <- table( cut( a, -3:2 ), b )
tabplot( tb )
tabplot( tb, rowlabs="right", col=heat.colors )

# Very similar plots
ptb <- sweep( tb, 2, apply( tb, 2, sum ), "/" )
par( mfrow=c(2,2) )
barplot( ptb, space=0 )
tabplot( tb, equal=TRUE, lwd=1 )
tabplot( tb, equal=TRUE, lwd=1, rowlabs="l" )
tabplot( tb, equal=FALSE, lwd=1, rowlabs="l" )
}
\keyword{hplot}
