\name{tbox}
\Rdversion{1.1}
\alias{tbox}
\alias{dbox}
\alias{fillarr}
\alias{boxarr}
\alias{boxes}
\alias{boxes.Lexis}
\alias{boxes.matrix}
\title{
  Draw boxes and arrows for illustration of multistate models.
  }
  \description{
  Boxes can be drawn with text (\code{tbox}) or a cross (\code{dbox}),
  and arrows pointing between the boxes (\code{boxarr}) can be drawn
  automatically not overlapping the boxes. The \code{boxes} method for
  \code{\link{Lexis}} objects generates displays of states with
  person-years and transitions with events or rates.
  }
\usage{
   tbox( txt, x, y, wd, ht,
         font=2, lwd=2,
         col.txt="black",
         col.border="black",
         col.bg="transparent" )
   dbox( x, y, wd, ht=wd,
         font=2, lwd=2, cwd=5,
         col.cross="black",
         col.border="black",
         col.bg="transparent"  )
   boxarr( b1, b2, offset=FALSE, pos=0.45, ... )
\method{boxes}{Lexis}( obj, file,
             detailed = FALSE,
               boxpos = FALSE,
                wmult = 1.5,
                hmult = 1.5*wmult,
                  cex = 1.5,
               show   = inherits( obj, "Lexis" ),
               show.Y = show,
              scale.Y = 1,
             digits.Y = 1,
               show.D = show,
              scale.D = FALSE,
             digits.D = as.numeric(as.logical(scale.D)),
                eq.wd = TRUE,
                eq.ht = TRUE,
                   wd,
                   ht,
               subset = NULL,
              exclude = NULL,
                 font = 2,
                  lwd = 2,
              col.txt = "black",
           col.border = col.txt,
               col.bg = "transparent",
              col.arr = "black",
              lwd.arr = 2,
             font.arr = 2,
              txt.arr = NULL,
          col.txt.arr = col.arr,
           offset.arr = 2, ... )
\method{boxes}{matrix}( obj, ... )
   fillarr( x1, y1, x2, y2, gap=2, fr=0.8,
            angle=17, lwd=2, length=par("pin")[1]/30, ... )
   }
\arguments{
  \item{txt}{Text to be placed inside the box.}
  \item{x}{x-coordinate of center of box.}
  \item{y}{y-coordinate of center of box.}
  \item{wd}{width of boxes in percentage of the plot width.}
  \item{ht}{height of boxes in percentage of the plot height.}
  \item{font}{Font for the text. Defaults to 2 (=bold).}
  \item{lwd}{Line width of the boxborders.}
  \item{col.txt}{Color for the text in boxes.}
  \item{col.border}{Color of the box border.}
  \item{col.bg}{Background color for the interior of the box.}
  \item{\dots}{Arguments to be passed on to the call of other functions.}
  \item{cwd}{Width of the lines in the cross.}
  \item{col.cross}{Color of the cross.}
  \item{b1}{Coordinates of the "from" box. A vector with 4 components,
            \code{x}, \code{y}, \code{w}, \code{h}.}
  \item{b2}{Coordinates of the "to" box; like \code{b1}.}
  \item{offset}{Logical. Should the arrow be offset a bit to the left.}
  \item{pos}{Numerical between 0 and 1, determines the position of the point
    on the arrow which is returned.}
  \item{obj}{A \code{\link{Lexis}} object or a transition matrix; that
    is a square matrix indexed by state in both dimensions, and the
    \eqn{(i,j)}th entry different from \code{NA} if a transition \eqn{i}
    to \eqn{j} can occur. If \code{show.D=TRUE}, the arrows between
    states are annotated by these numbers. If \code{show.Y=TRUE}, the
    boxes representing states are annotated by the numbers in the
    diagonal of \code{obj}.}
  \item{file}{Name of the file with the code reproducing the plot.}
  \item{detailed}{Should the output of R-code be detailed, showing all parameters?}
  \item{boxpos}{If \code{TRUE} the boxes are positioned  equidistantly on a
                circle, if \code{FALSE} (the default) you are queried to
                click on the screen for the positions. This argument can also
                be a named list with elements \code{x} and
                \code{y}, both numerical vectors, giving the centers of
                the boxes.}
  \item{wmult}{Multiplier for the width of the box relative to the width of the
               text in the box.}
  \item{hmult}{Multiplier for the height of the box relative to the height of the
               text in the box.}
  \item{cex}{Character expansion for text in the box.}
  \item{show}{Should person-years and transitions be put in the plot.
              Ignored if \code{obj} is not a \code{Lexis} object.}
  \item{show.Y}{If logical: Should person-years be put in the boxes.
              If numeric: Numbers to put in boxes.}
  \item{scale.Y}{What scale should be used for annotation of person-years.}
  \item{digits.Y}{How many digits after the decimal point should be used for the
                  person-years.}
  \item{show.D}{Should transitions be put alongside the arrows.
              Ignored if \code{obj} is not a \code{Lexis} object.}
  \item{scale.D}{If this a scalar, rates instead of no. transitions are printed
                 at the arrows, scaled by \code{scale.D}.}
  \item{digits.D}{How many digits after the decimal point should be used for the
                  rates.}
  \item{eq.wd}{Should boxes all have the same width?}
  \item{eq.ht}{Should boxes all have the same height?}
  \item{subset}{Draw only boxes and arrows for a subset of the states.
                Can be given either as a numerical vector or character
                vector state names.}
  \item{exclude}{Exclude states from the plot. The complementary of \code{subset}.
    Ignored if \code{subset} is given.}
  \item{col.arr}{Color of the arrows between boxes.
    A vector of character strings, the arrows are referred to as the
    row-wise sequence of non-NA elements of the transition matrix.
    Thus the first ones refer to the transitions out of state 1, in
    order of states.}
  \item{lwd.arr}{Line withs of the arrows.}
  \item{font.arr}{Font of the text annotation the arrows.}
  \item{txt.arr}{Text put on the arrows.}
  \item{col.txt.arr}{Colors for text on the arrows.}
  \item{offset.arr}{The amount offset between arrows representing
    two-way transitions, that is where there are arrows both ways
    between two boxes.}
  \item{x1}{x-coordinate of the starting point.}
  \item{y1}{y-coordinate of the starting point.}
  \item{x2}{x-coordinate of the end point.}
  \item{y2}{y-coordinate of the end point.}
  \item{gap}{Length of the gap between the box and the ends of the arrows.}
  \item{fr}{Length of the arrow as the fraction of the distance between the
            boxes. Ignored unless given explicitly, in which case any value
            given for \code{gap} is ignored.}
  \item{angle}{What angle should the arrow-head have?}
  \item{length}{Length of the arrow head in inches. Defaults to 1/30 of the
                physical width of the plot.}
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
  area (with axes 0 to 100 in both directions) and if
  \code{boxpos=FALSE} (the default) prompts you to click on the
  locations for the state boxes, and then draws arrows implied by the
  actual transitions in the \code{Lexis} object.

  A transition matrix can also be supplied, in which case the row/column
  names are used as state names, diagnonal elements taken as
  person-years, and off-diagnonal elements as number of transitions.
  This also works for \code{boxes.matrix}.

  Optionally returns the R-code reproducing the plot in a file, which
  can be useful if you want to produce exactly the same plot with
  differing arrow colors etc.

  \code{boxarr} draws an arrow between two boxes, on the line connecting
  the two box centers. The \code{offset} argument is used to offset the
  arrow a bit to the left (as seen in the direction of the arrow) on order
  to accommodate arrows both ways between boxes. \code{boxarr} returns a named
  list with elements \code{x}, \code{y} and \code{d}, where the two former
  give the location of a point on the arrow used for printing (see argument
  \code{pos}) and the latter is a unit vector in the
  direction of the arrow, which is used by \code{boxes.Lexis} to
  position the annotation of arrows with the number of transitions.
  \code{fill.arr} is just a utility drawing nicer arrows than the default
  \code{\link{arrows}} command, basically by using filled arrow-heads; called
  by \code{boxarr}.
  }
\value{The functions \code{tbox} and \code{dbox} return the location and
       dimension of the boxes, \code{c(x,y,w,h)}, which are designed to be used
       as input to the \code{boxarr} function.

       The \code{boxarr} function returns the coordinates (as a named
       list with names \code{x} and \code{y}) of a point on the
       arrow, designated to be used for annotation of the arrow.

       The function \code{boxes.Lexis} returns (invisibly) a list with
       two elements: \code{Boxes} - a dataframe with one row
       per box and columns
       \code{pl.nam},
       \code{xx},
       \code{yy},
       \code{wd},
       \code{ht},
       \code{font},
       \code{lwd},
       \code{col.txt},
       \code{col.border} and
       \code{col.bg},
       and \code{Tmat}, the transition matrix.
       }
\author{Bendix Carstensen}
\examples{
par( mar=c(0,0,0,0), cex=1.5 )
plot( NA,
      bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
bw  <- tbox( "Well"    , 10, 60, 22, 10, col.txt="blue" )
bo  <- tbox( "other Ca", 45, 80, 22, 10, col.txt="gray" )
bc  <- tbox( "Ca"      , 45, 60, 22, 10, col.txt="red" )
bd  <- tbox( "DM"      , 45, 40, 22, 10, col.txt="blue" )
bcd <- tbox( "Ca + DM" , 80, 60, 22, 10, col.txt="gray" )
bdc <- tbox( "DM + Ca" , 80, 40, 22, 10, col.txt="red" )
      boxarr( bw, bo , col=gray(0.7), lwd=3 )
# Note the argument adj= can takes values outside (0,1)
text( boxarr( bw, bc , col="blue", lwd=3 ),
      expression( lambda[Well] ), col="blue", adj=c(1,-0.2), cex=0.8 )
      boxarr( bw, bd , col=gray(0.7) , lwd=3 )
      boxarr( bc, bcd, col=gray(0.7) , lwd=3 )
text( boxarr( bd, bdc, col="blue", lwd=3 ),
      expression( lambda[DM] ), col="blue", adj=c(1.1,-0.2), cex=0.8 )

# Set up a transition matrix allowing recovery
tm <- rbind( c(NA,1,1), c(1,NA,1), c(NA,NA,NA) )
rownames(tm) <- colnames(tm) <- c("Cancer","Recurrence","Dead")
boxes.Lexis( tm, file="", boxpos=TRUE )
boxes.Lexis( tm, file="", detailed=TRUE, boxpos=TRUE )
# Illustrate texting of arrows
boxes.Lexis( tm, boxpos=TRUE, txt.arr=c("en","to","tre","fire") )
boxes.Lexis( tm, boxpos=TRUE, txt.arr=c(expression(lambda[C]),
                                        expression(mu[C]),
                                        "recovery",
                                        expression(mu[R]) ) )

# Set up a Lexis object
data(DMlate)
str(DMlate)
dml <- Lexis( entry=list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit=list(Per=dox),
        exit.status=factor(!is.na(dodth),labels=c("DM","Dead")),
               data=DMlate[1:1000,] )
# Split follow-up at Insulin
dmi <- cutLexis( dml, cut=dml$doins, new.state="Ins", pre="DM" )
summary( dmi )
boxes( dmi, boxpos=TRUE )
# Set up a bogus recovery date
dmi$dorec <- dmi$doins + runif(nrow(dmi),0.5,10)
dmi$dorec[dmi$dorec>dmi$dox] <- NA
dmR <- cutLexis( dmi, cut=dmi$dorec, new.state="DM", pre="Ins" )
summary( dmR )
boxes( dmR, boxpos=TRUE )
    }
\seealso{
   \code{\link{tmat.Lexis}}
   }
\keyword{survival}
\keyword{hplot}
\keyword{iplot}
