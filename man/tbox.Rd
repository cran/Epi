\name{tbox}
\Rdversion{1.1}
\alias{tbox}
\alias{dbox}
\alias{fillarr}
\alias{boxarr}
\alias{boxes.Lexis}
\title{
  Draw boxes and arrows for illustration of multistate models.
  }
\description{
  Boxes can be drawn with text (\code{tbox}) or a cross (\code{dbox}), and
  arrows pointing between the boxes can be drawn automatically not overlapping
  the \code{boxarr} boxes.
  }
\usage{
   tbox( txt, x, y, w, h,
         font = 2, txt.col = "black", lwd = 2,
         border = "black", col = "transparent")
   dbox( x, y, w, h=w,
         font=2, cross.col="black", cwd=5,
         lwd=2, border="black", col="transparent" )
   boxarr( b1, b2, ... )
   \method{boxes}{Lexis}( obj, file, mult=1.5, cex=1.5 )
   fillarr( x1, y1, x2, y2, fr=0.8,
            angle=17, lwd=2, length=par("pin")[1]/30, ... )
   }
\arguments{
  \item{txt}{Text to be placed inside the box.}
  \item{x}{x-coordinate of center of box.}
  \item{y}{y-coordinate of center of box.}
  \item{w}{width of box.}
  \item{h}{height of box.}
  \item{font}{Font for the text.}
  \item{txt.col}{Color for the text.}
  \item{lwd}{Line width of the box / arrow.}
  \item{border}{Color of the box border.}
  \item{col}{Background color for the interior of the box.}
  \item{cross.col}{Color of the cross.}
  \item{cwd}{Width of the lines in the cross.}
  \item{b1}{Coordinates of the "from" box. A vector with 4 components,
    \code{x}, \code{y}, \code{w}, \code{h}.}
  \item{b2}{Coordinates of the "to" box.}
  \item{...}{Arguments to be passed on to the call of \code{\link{arrows}}.} 
  \item{obj}{A \code{\link{Lexis}} object, or a transition matrix; that
    is a matrix}
  \item{file}{Name of the file with the code reproducing the plot.}
  \item{mult}{Multiplier for the withdt of the box relative to the width of the
              text in the box.}
  \item{cex}{Character expansion for text in the box.}
  \item{x1}{x-coordinate of the starting point.}
  \item{y1}{y-coordinate of the starting point.}
  \item{x2}{x-coordinate of the ending point.}
  \item{y2}{y-coordinate of the ending point.}
  \item{fr}{Length of the arrow as the fraction of the distance between the
            boxes.}
  \item{angle}{What angle should the arrow-head have?}
  \item{length}{Length of the arrow head in inches.}
  }
  \details{    
  These functions are designed to facilitate the drawing of multistate
  models, mainly by automatic calculation of the arrows between boxes.
  \code{tbox} draws a box with centered text, and returns a vector of
  location, height and width of the box. This is used when drawing
  arrows between boxes. \code{dbox} draws a box with a cross,
  symbolizing a death state. \code{boxarr} draws an arrow between two
  boxes, making sure it does not intersect the boxes. Only straight
  lines are drawn. 
  
  \code{boxes.Lexis} takes as input a Lexis object sets up an empty plot
  area and prompts you to click on the locations for the state boxes, and the
  draws the model. A transition matrix can also be supplied, in which
  case the row/column names are used as state names. Optionally returns
  a file with R-code reproducing the plot.             
  }
\value{The functions \code{tbox} and \code{dbox} return the location and
       dimension of the boxes, \code{c(x,y,w,h)}, which are designed to be used
       as input to the \code{boxarr} function.

       The \code{boxarr} function returns the coordinates (as a named
       list with names \code{x} and \code{y}) of the midpoint of the
       arrow, designated to be used for annotation of the arrow.
       } 
\author{Bendix Carstensen}
\examples{
par( mar=c(0,0,0,0), cex=1.5 )
plot( NA,
      bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
bw  <- tbox( "Well"    , 10, 60, 22, 10, txt.col="blue" )
bo  <- tbox( "other Ca", 45, 80, 22, 10, txt.col="gray" )
bc  <- tbox( "Ca"      , 45, 60, 22, 10, txt.col="red" )
bd  <- tbox( "DM"      , 45, 40, 22, 10, txt.col="blue" )
bcd <- tbox( "Ca + DM" , 80, 60, 22, 10, txt.col="gray" )
bdc <- tbox( "DM + Ca" , 80, 40, 22, 10, txt.col="red" )
      boxarr( bw, bo , col=gray(0.7), lwd=3 )
# Note the argument adj= can takes values outside (0,1) 
text( boxarr( bw, bc , col="blue", lwd=3 ),
      expression( lambda[Well] ), col="blue", adj=c(1,-0.2), cex=0.8 )
      boxarr( bw, bd , col=gray(0.7) , lwd=3 )
      boxarr( bc, bcd, col=gray(0.7) , lwd=3 )
text( boxarr( bd, bdc, col="blue", lwd=3 ),
      expression( lambda[DM] ), col="blue", adj=c(1.1,-0.2), cex=0.8 )
# Set up a transition matrix
tm <- rbind( c(NA,1,1), c(NA,NA,1), c(NA,NA,NA) )
rownames(tm) <- colnames(tm) <- c("Cancer","Recurrence","Dead")
\dontrun{boxes.Lexis( tm )}
    }
\keyword{survival}
\keyword{hplot}
\keyword{iplot}
