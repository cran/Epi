tbox <-
function( txt, x, y, w, h,
          font=2, txt.col="black",
          lwd=2, border="black", col="transparent", ... )
{
rect( x-w/2, y-h/2, x+w/2, y+h/2, lwd=lwd, border=border, col=col )
text( x, y, txt, font=font, col=txt.col, ... )
invisible( c( x, y, w, h ) )
}

dbox <-
function( x, y, w, h=w,
          font=2, cross.col="black", cwd=5,
          lwd=2, border="black", col="transparent" )
{
rect( x-w/2, y-h/2, x+w/2, y+h/2, lwd=lwd, border=border, col=col )
ch <- h*2/3
segments( c(x     , x-ch/3),
          c(y+ch/2, y+ch/6),
          c(x     , x+ch/3),
          c(y-ch/2, y+ch/6), lwd=cwd, col=cross.col )
invisible( c( x, y, w, h ) )
}

fillarr <-
function( x1, y1, x2, y2, gap=2, fr=0.8,
          angle=17, lwd=2, length=par("pin")[1]/30, ... )
{
fr <- 1-gap/sqrt((x1-x2)^2+(y1-y2)^2)
if( !missing(fr)  ) if( fr > 1 ) fr <- fr/100
for( a in 1:angle )
arrows( x1 + (x2-x1)*(1-fr)/2,
        y1 + (y2-y1)*(1-fr)/2,
        x2 - (x2-x1)*(1-fr)/2,
        y2 - (y2-y1)*(1-fr)/2, angle=a, lwd=lwd, ... )
}

std.vec <-
function( a, b )
{
  l <- sqrt(a^2+b^2)
  if( l==0 )
  return( c(0,0) )
  else
  return( c(a/l,b/l) )
}

boxarr <-
function( b1, b2, offset=FALSE, pos=0.6, ... )
{
# If we want to offset the arrow a bit to the left
if( offset ) d <- std.vec( b2[1]-b1[1], b2[2]-b1[2] )
else d <- c(0,0)
x1 <- b1[1] - d[2]
y1 <- b1[2] + d[1]
w1 <- b1[3]
h1 <- b1[4]
x2 <- b2[1] - d[2]
y2 <- b2[2] + d[1]
w2 <- b2[3]
h2 <- b2[4]
hx1 <- x1 + ifelse( (y2-y1)!=0, (x2-x1)*((h1/2)/abs(y2-y1)), sign(x2-x1)*w1/2 )
vx1 <- x1 + ifelse( (x2-x1)!=0, (x2-x1)*((w1/2)/abs(x2-x1)), 0    )
hx2 <- x2 + ifelse( (y1-y2)!=0, (x1-x2)*((h2/2)/abs(y1-y2)), sign(x1-x2)*w2/2 )
vx2 <- x2 + ifelse( (x1-x2)!=0, (x1-x2)*((w2/2)/abs(x1-x2)), 0    )
hy1 <- y1 + ifelse( (y2-y1)!=0, (y2-y1)*((h1/2)/abs(y2-y1)), 0    )
vy1 <- y1 + ifelse( (x2-x1)!=0, (y2-y1)*((w1/2)/abs(x2-x1)), sign(y2-y1)*h1/2 )
hy2 <- y2 + ifelse( (y1-y2)!=0, (y1-y2)*((h2/2)/abs(y1-y2)), 0    )
vy2 <- y2 + ifelse( (x1-x2)!=0, (y1-y2)*((w2/2)/abs(x1-x2)), sign(y1-y2)*h2/2 )
if( abs(vy1-y1) < h1/2 ) { bx1 <- vx1
                           by1 <- vy1 }
                    else { bx1 <- hx1
                           by1 <- hy1 }
if( abs(vy2-y2) < h2/2 ) { bx2 <- vx2
                           by2 <- vy2 }
                    else { bx2 <- hx2
                           by2 <- hy2 }
fillarr( bx1, by1, bx2, by2, ... )
invisible( list( x = bx1*(1-pos)+bx2*pos,
                 y = by1*(1-pos)+by2*pos,
                 d = std.vec(bx2-bx1,by2-by1) ) )
}

boxes <- function (obj, ...) UseMethod("boxes")

boxes.Lexis <-
function( obj, file,
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
                eq.ht = TRUE, ... )
{
if( inherits(obj,"Lexis") ) tm <- tmat( obj )
else if( is.matrix(obj) & diff(dim(obj))==0 ) tm <- obj
else stop( "First argument must be a Lexis object or a square matrix.\n" )

# First get the number of transitions, then write it all to a file
# and then source it to do the plot
                      st.nam <- colnames( tm )
if( is.null(st.nam) ) st.nam <- paste(1:ncol(tm))
            pl.nam <- st.nam
      n.st <- length( st.nam )

# Do we want to show person-years and events
if( show & inherits( obj, "Lexis" ) )
  {
  SM <- summary(obj,simplify=FALSE,scale=scale.Y)
  Y <- SM[[1]][-n.st-1,"Risk time:"]
  D <- SM[[1+as.logical(scale.D)]][1:n.st,1:n.st] * ifelse(scale.D,scale.D,1)
  }
# No extra line with person-years when they are NA
if( show.Y ) pl.nam <- gsub("\\\nNA", "",
    paste( st.nam, formatC( Y, format="f", digits=digits.Y, big.mark="," ), sep="\n" ) )

# Her comes the plot
par( mar=c(0,0,0,0), cex=cex )
plot( NA,
      bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
# String height and width only meaningful after a plot has been called
ht <- strheight( pl.nam ) * hmult
wd <- strwidth(  pl.nam ) * wmult
if( eq.ht ) ht <- rep( max(ht), length(ht) )
if( eq.wd ) wd <- rep( max(wd), length(wd) )

# If not supplied, ask for positions of boxes
if( is.list(boxpos) )
  {
  if( names(boxpos) != c("x","y") )
    stop( "The list given in 'boxpos=' must have names 'x' and 'y'" )
  if( length(boxpos$x) != n.st | length(boxpos$y) != n.st )
    stop( "The elements of boxpos must both have length equal to no. states", n.st )
  xx <- boxpos$x
  yy <- boxpos$y
  }
if( is.logical(boxpos) )
  if( boxpos )
  {
  ang <- pi - 2*pi*((1:n.st-0.5)/n.st)
  xx <- cos( ang ) * 35 + 50
  yy <- sin( ang ) * 35 + 50
  }
  else
  {
  xx <- yy <- numeric(n.st)
  for( i in 1:n.st )
     {
     cat( "\nClick for level ", st.nam[i] )
     flush.console()
     pt <- locator(1)
     xx[i] <- pt$x
     yy[i] <- pt$y
     tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i] )
     }
  cat( "\n" )
  }
# Plot the boxes and record position and size
b <- list()
for( i in 1:n.st ) b[[i]] <- tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i] )
for( i in 1:n.st ) for( j in 1:n.st )
  {
  if( !is.na(tm[i,j]) & i!=j )
    {
    arr <- boxarr( b[[i]], b[[j]], offset=!is.na(tm[j,i]), ... )
    if( show.D )
    text( arr$x-arr$d[2], arr$y+arr$d[1],
          formatC( D[i,j], format="f", digits=digits.D ),
          adj=as.numeric(c(arr$d[2]>0,arr$d[1]<0)),
          font=2, col="black" )
    }
  }
# Redraw the boxes with white background
for( i in 1:n.st ) tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i], col="white" )

#################################################################
# If we want the code we just print the entire function code,
# where we have replaced the query for the coordinates of the
# boxes with the values we obtained from locator()
if( !missing(file) )
{
xx <- round(xx)
yy <- round(yy)

# Here we start writing the file, beginning with setting the parameters
cat( '
obj <- ', deparse( substitute( obj ) ), '\n
if( inherits(obj,"Lexis") ) tm <- tmat.Lexis( obj ) else tm <- obj
### Position of the boxes:
      xx <- c(', paste( xx, collapse=", " ),')
      yy <- c(', paste( yy, collapse=", " ),')
     cex <-', cex, '      # How should text and numbers be scaled
   wmult <-', wmult, '    # Extra box-width relative to string width
   hmult <-', hmult, '    # Extra box-height relative to string height
   eq.wd <-', eq.wd, '    # All boxes the same width
   eq.ht <-', eq.ht, '    # All boxes the same height
  show.Y <-', show.Y, '   # Show number of person-years in boxes
 scale.Y <-', scale.Y, '  # How should person-years be scaled
digits.Y <-', digits.Y, ' # How should person-years be printed
  show.D <-', show.D, '   # Show number of events on arrows
 scale.D <-', scale.D, '  # How should rates be scaled
digits.D <-', digits.D, ' # How should rates be printed
                      st.nam <- colnames( tm )
if( is.null(st.nam) ) st.nam <- paste(1:ncol(tm))
            pl.nam <- st.nam
      n.st <- length( st.nam )

# If we want to show person-years and events / rates
  SM <- summary(obj,simplify=FALSE,scale=scale.Y)
  Y <- SM[[1]][-n.st-1,"Risk time:"]
  D <- SM[[1+as.logical(scale.D)]][1:n.st,1:n.st] * ifelse(scale.D,scale.D,1)
# No extra line with person-years when they are NA
# and adjust the printing of the person-Years
if( show.Y )
  {
  pl.nam <- gsub( "\\\nNA", "",
                  paste( st.nam,
                         formatC( Y, format="f", digits=',digits.Y,',
                                  big.mark=","),
                         sep="\n" ) )
  }
  
#################################################################
# Here comes the plot
par( mar=c(0,0,0,0), cex=cex )
plot( NA, bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
# String height and width only meaningful after a plot has been called
ht <- strheight( pl.nam ) * hmult
wd <- strwidth(  pl.nam ) * wmult

# Should all boxes be of same height / width:
if( eq.ht ) ht <- rep( max(ht), length(ht) )
if( eq.wd ) wd <- rep( max(wd), length(wd) )

# Plot the boxes
b <- list()
for( i in 1:n.st )
   b[[i]] <- tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i] )

# Plot the arrows and the text along them
for( i in 1:n.st ) for( j in 1:n.st )
   {
   if( !is.na(tm[i,j]) )
     {
     arr <- boxarr( b[[i]], b[[j]], offset=!is.na(tm[j,i]) )
     if( show.D ) text( arr$x-arr$d[2], arr$y+arr$d[1],
                        formatC( D[i,j], format="f", digits=digits.D ),
                        adj=as.numeric(c(arr$d[2]>0,arr$d[1]<0)),
                        font=2, col="black" )
     }
   }
# Redraw the boxes with white background to move the arrows behind them
for( i in 1:n.st ) tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i], col="white" )
     ', file=file )
}
}
