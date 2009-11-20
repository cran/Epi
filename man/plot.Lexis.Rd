\name{plot.Lexis}
\alias{plot.Lexis}
\alias{points.Lexis}
\alias{lines.Lexis}
\alias{PY.ann}
\alias{PY.ann.Lexis}
\title{Lexis diagrams}
\description{
  The follow-up histories represented by a Lexis object can be plotted
  using one or two dimensions.  The two dimensional plot is a Lexis
  diagram showing follow-up time simultaneously on two time scales.
}
\usage{
\method{plot}{Lexis}(x, time.scale = NULL, type="l", breaks="lightgray", ...)
\method{points}{Lexis}(x, time.scale = options()[["Lexis.time.scale"]] , ...)
\method{lines}{Lexis}(x, time.scale = options()[["Lexis.time.scale"]], ...)
\method{PY.ann}{Lexis}(x, time.scale = options()[["Lexis.time.scale"]], digits=1, ...)
}
\arguments{
  \item{x}{An object of class \code{Lexis}}
  \item{time.scale}{A vector of length 1 or 2 giving the time scales to
    be plotted either by name or numerical order}
  \item{type}{Character indication what to draw: "n" nothing (just set up the
    diagram), "l" - liefelines, "p" - endpoints of follow-up, "b" - both
    lifelines and endpoints.}
  \item{breaks}{a string giving the colour of grid lines to be drawn
    when plotting a split Lexis object. Grid lines can be suppressed by
    supplying the value \code{NULL} to the \code{breaks} argument}
  \item{digits}{Numerical. How many digits after the demimal points should be
       when plotting the person-years.}
  \item{\dots}{Further graphical parameters to be passed to the plotting
    methods.

    Grids can be drawn (behind the life lines) using the following
    parameters in \code{plot}: 
    \itemize{
      \item \code{grid} If logical, a background grid is set up
      using the axis ticks. If a list, the first component is used as
      positions for the vertical lines and the last as positions for the
      horizontal. If a nunerical vector, grids on both axes are set up
      using the distance between the numbers.
      \item \code{col.grid="lightgray"} Color of the background grid.
      \item \code{lty.grid=2} Line type for the grid.
      \item \code{coh.grid=FALSE} Should a 45 degree grid be plotted?}
} }
\details{
  The plot method for \code{Lexis} objects traces ``life lines'' from
  the start to the end of follow-up.  The \code{points} method plots
  points at the end of the life lines.

  If \code{time.scale} is of length 1, the life lines are drawn
  horizontally, with the time scale on the X axis and the id value on the Y
  axis. If \code{time.scale} is of length 2, a Lexis diagram is
  produced, with diagonal life lines plotted against both time scales
  simultaneously.

  If \code{lex} has been split along one of the time axes by a call to
  \code{splitLexis}, then vertical or horizontal grid lines are plotted
  (on top of the life lines) at the break points.
  
  \code{PY.ann} writes the length of each (segment of) life line at the middle
  of the line. Not advisable to use with large cohorts. Another example is
  in the example file for \code{\link{occup}}.
}
\author{Martyn Plummer}
\examples{
# A small bogus cohort
xcoh <- structure( list( id = c("A", "B", "C"),
                      birth = c("14/07/1952", "01/04/1957", "10/06/1987"),
                      entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                       exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                       fail = c(1, 0, 1) ),
                     .Names = c("id", "birth", "entry", "exit", "fail"),
                  row.names = c("1", "2", "3"),
                      class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh$bt <- cal.yr( xcoh$birth, format="\%d/\%m/\%Y" )
xcoh$en <- cal.yr( xcoh$entry, format="\%d/\%m/\%Y" )
xcoh$ex <- cal.yr( xcoh$exit , format="\%d/\%m/\%Y" )

# See how it looks
xcoh

# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( per=en ),
                exit = list( per=ex, age=ex-bt ),
         exit.status = fail,
                data = xcoh )

# Default plot of follow-up
plot( Lcoh )
# Show follow-up time
PY.ann( Lcoh )

# Show exit status
plot( Lcoh, type="b" )
# Same but failures only
plot( Lcoh, type="b", pch=c(NA,16)[Lcoh$fail+1] )

# With a grid and deaths as endpoints
plot( Lcoh, grid=0:10*10, col="black" )
points( Lcoh, pch=c(NA,16)[Lcoh$lex.Xst+1] )
# With a lot of bells and whistles:
plot( Lcoh, grid=0:20*5, col="black", xaxs="i", yaxs="i",
      xlim=c(1960,2010), ylim=c(0,50), lwd=3, las=1 )
points( Lcoh, pch=c(NA,16)[Lcoh$lex.Xst+1], col="red", cex=1.5 )
}
\seealso{\code{\link{Lexis}}, \code{\link{splitLexis}}}
\keyword{hplot}
\keyword{aplot}
