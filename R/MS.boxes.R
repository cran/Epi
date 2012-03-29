tbox <-
function( txt, x, y, wd, ht,
          font=2, lwd=2,
          col.txt="black",
          col.border="black",
          col.bg="transparent" )
{
rect( x-wd/2, y-ht/2, x+wd/2, y+ht/2,
      lwd=lwd, border=col.border, col=col.bg )
text( x, y, txt, font=font, col=col.txt )
invisible( c( x, y, wd, ht ) )
}

dbox <-
function( x, y, wd, ht=wd,
          font=2, lwd=2, cwd=5,
          col.cross="black",
          col.border="black",
          col.bg="transparent" )
{
rect( x-wd/2, y-ht/2, x+wd/2, y+ht/2, lwd=lwd, border=col.border, col=col.bg )
ch <- ht*2/3
segments( c(x     , x-ch/3),
          c(y+ch/2, y+ch/6),
          c(x     , x+ch/3),
          c(y-ch/2, y+ch/6), lwd=cwd, col=col.cross )
invisible( c( x, y, wd, ht ) )
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
function( b1, b2, offset=FALSE, pos=0.45, ... )
{
# If we want to offset the arrow a bit to the left, we compute
# A unit vector in the direction of the arrow and add twice
# the orthogonal of this to the coordinates
d  <- std.vec( b2[1]-b1[1], b2[2]-b1[2] )
dd <- d * offset
x1 <- b1[1] - dd[2]
y1 <- b1[2] + dd[1]
w1 <- b1[3]
h1 <- b1[4]
x2 <- b2[1] - dd[2]
y2 <- b2[2] + dd[1]
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
                 d = d ) )
}

boxes <- function (obj, ...) UseMethod("boxes")

boxes.matrix <-
function( obj, ... )
  {
  Epi:::boxes.Lexis( obj, ... )
  }

boxes.Lexis <-
function( obj, file, detailed=FALSE,
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
{
if( inherits(obj,"Lexis") )
  {
  if( !is.factor(obj$lex.Cst) | !is.factor(obj$lex.Xst) ) obj <- factorize( obj )
  tm <- tmat( obj, Y=TRUE )
  }
else if( is.matrix(obj) & diff(dim(obj))==0 )
  {
  tm <- obj
  }
else stop( "First argument must be a Lexis object or a square matrix.\n" )

# Put the transitions into D and the diagnonal into Y.
D <- tm
diag( D ) <- NA
Y <- diag( tm )

# Derive state names, no. states and no. transitions
                      st.nam <- colnames( tm )
if( is.null(st.nam) ) st.nam <- paste(1:ncol(tm))
            pl.nam <- st.nam
      n.st <- length( st.nam )
      n.tr <- sum( !is.na(tm) )

# If we want to show person-years and events we compute them
if( inherits(obj,"Lexis") )
  {
  if( show )
    {
    SM <- summary(obj,simplify=FALSE,scale=scale.Y)
    Y <- SM[[1]][1:n.st,"Risk time:"]
    D <- SM[[1+as.logical(scale.D)]][1:n.st,1:n.st] * ifelse(scale.D,scale.D,1)
    }
  }

# Explicitly given numbers in boxes ?
if( is.numeric(show.Y) )
  {
  Y <- show.Y
  show.Y <- TRUE
  }

# No extra line with person-years when they are NA
if( show.Y ) pl.nam <- gsub( "\\\nNA",
                             "",
                             paste( st.nam,
                                    formatC( Y,
                                             format="f",
                                             digits=digits.Y,
                                             big.mark="," ),
                                    sep="\n" ) )

# Any subsetting:
sbst <- 1:nrow(tm)
if( !is.null(exclude) )
  {
  if( is.character(exclude) )
    exclude <- match( exclude, rownames(tm) )
  sbst <- sbst[-exclude]
  }
if( !is.null(subset) )
  {
  if( is.character(subset) )
      sbst <- match( subset, rownames(tm) )
  else sbst <- subset
  }
subset <- sbst

# Recycling of box-arguments
if( !missing(ht) )
if( length(ht         )<n.st ) ht         <- rep(ht        ,n.st)[1:n.st]
if( !missing(wd) )
if( length(wd         )<n.st ) wd         <- rep(wd        ,n.st)[1:n.st]
if( length(font       )<n.st ) font       <- rep(font      ,n.st)[1:n.st]
if( length(lwd        )<n.st ) lwd        <- rep(lwd       ,n.st)[1:n.st]
if( length(col.border )<n.st ) col.border <- rep(col.border,n.st)[1:n.st]
if( length(col.bg     )<n.st ) col.bg     <- rep(col.bg    ,n.st)[1:n.st]
if( length(col.txt    )<n.st ) col.txt    <- rep(col.txt   ,n.st)[1:n.st]
# Recycling of arrow-arguments
if( length(col.arr    )<n.tr ) col.arr    <- rep(col.arr    ,n.tr)[1:n.tr]
if( length(col.txt.arr)<n.tr ) col.txt.arr<- rep(col.txt.arr,n.tr)[1:n.tr]
if( length(lwd.arr    )<n.tr ) lwd.arr    <- rep(lwd.arr    ,n.tr)[1:n.tr]
if( length(font.arr   )<n.tr ) font.arr   <- rep(font.arr   ,n.tr)[1:n.tr]

# Here comes the plot
# First setting up the plot area, and restoring the plot parameters later
opar <- par( mar=c(0,0,0,0), cex=cex )
on.exit( par(opar) )

plot( NA,
      bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
# String height and width only meaningful after a plot has been called
if( missing(ht) )
  {
  ht <- strheight( pl.nam ) * hmult
  if( eq.ht ) ht <- rep( max(ht), length(ht) )
  }
if( missing(wd) )
  {
  wd <- strwidth(  pl.nam ) * wmult
  if( eq.wd ) wd <- rep( max(wd), length(wd) )
  }

# If not supplied, ask for positions of boxes
if( is.list(boxpos) )
  {
  names(boxpos) <- tolower( names(boxpos) )
  if( length(intersect(names(boxpos),c("x","y")))<2 )
    stop( "The list given in 'boxpos=' must have components 'x' and 'y'" )
  if( length(boxpos$x) != n.st | length(boxpos$y) != n.st )
    stop( "The elements 'x' and 'y' of boxpos must both have length equal to no. states", n.st )
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
  for( i in subset )
     {
     cat( "\nClick for level ", st.nam[i] )
     flush.console()
     pt <- locator(1)
     xx[i] <- pt$x
     yy[i] <- pt$y
     tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i], ... )
     }
  cat( "\n" )
  }
# Plot the boxes and record position and size
b <- list()
for( i in subset ) b[[i]] <- tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i],
                                font=font[i],
                                  lwd=lwd[i],
                          col.txt=col.txt[i],
                    col.border=col.border[i],
                            col.bg=col.bg[i] )
# Arrows and text on them
for( i in subset ) for( j in subset )
  {
  if( !is.na(tm[i,j]) & i!=j )
    {
    arrowtext <- NULL
    # Which number of arrow is currently processed?
    a <- sum(!is.na(tm[1:i,])) - sum(!is.na(tm[i,j:n.st])) + 1
    arr <- boxarr( b[[i]], b[[j]],
                   offset=(!is.na(tm[j,i]))*offset.arr,
                   lwd=lwd.arr[a], col=col.arr[a], ... )
    if( show.D & D[i,j]>0 )
      arrowtext <- formatC( D[i,j], format="f", digits=digits.D, big.mark="," )
    else
    if( !is.null(txt.arr) )
      arrowtext <- txt.arr[a]
    if( !is.null(arrowtext) )
    text( arr$x-arr$d[2], arr$y+arr$d[1],
          arrowtext,
          adj=as.numeric(c(arr$d[2]>0,arr$d[1]<0)),
          font=font.arr[a], col=col.txt.arr[a] )
    }
  }
# Redraw the boxes with white background to remove any arrows
for( i in subset ) tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i],
                         lwd=lwd[i], col.bg="white" )
# Then redraw the boxes again
for( i in subset ) tbox( pl.nam[i], xx[i], yy[i], wd[i], ht[i],
                                font=font[i],
                                  lwd=lwd[i],
                          col.txt=col.txt[i],
                    col.border=col.border[i],
                            col.bg=col.bg[i] )


#################################################################
# If we want the code we just print the entire function code,
# where we have replaced the query for the coordinates of the
# boxes with the values we obtained from locator()
if( !missing(file) )
{
xx <- round(xx)
yy <- round(yy)
# Here we write the code-file which will reproduce the plot

# Simple plot
if( !detailed )
cat('
boxes.Lexis( ', deparse( substitute( obj ) ),',
           boxpos = list( x=c(', paste( xx, collapse=", " ),'),
                          y=c(', paste( yy, collapse=", " ),') ),
              cex =', cex,',\t # How should text and numbers be scaled
            wmult =', wmult,',\t # Extra box-width relative to string width
            hmult =', hmult,',\t # Extra box-height relative to string height
            eq.wd =', eq.wd,',\t # All boxes the same width
            eq.ht =', eq.ht,',\t # All boxes the same height
           show.Y =', show.Y,',\t # Show number of person-years in boxes
          scale.Y =', scale.Y,',\t\t # How should person-years be scaled
         digits.Y =', digits.Y,',\t\t # How should person-years be printed
           show.D =', show.D,',\t # Show number of events on arrows
          scale.D =', scale.D,',\t # How should rates be scaled
         digits.D =', digits.D,')\t\t # How should rates be printed \n',
 file=file )
else
cat('
boxes.Lexis( ', deparse( substitute( obj ) ),',
           boxpos = list( x=c(', paste( xx, collapse="," ),'),
                          y=c(', paste( yy, collapse="," ),') ),
              cex =', cex,',\t # How should text and numbers be scaled
               wd = c(', paste( round(wd), collapse="," ),'),\t # Width of boxes
               ht = c(', paste( round(ht), collapse="," ),'),\t # Height of boxes
           show.Y =', show.Y,',\t # Show number of person-years in boxes
          scale.Y =', scale.Y,',\t\t # How should person-years be scaled
         digits.Y =', digits.Y,',\t\t # How should person-years be printed
           show.D =', show.D,',\t # Show number of events on arrows
          scale.D =', scale.D,',\t # How should rates be scaled
         digits.D =', digits.D,',\t\t # How should rates be printed
          col.txt = c("', paste(col.txt   ,collapse="\",\""), '"),\t # color of the text
       col.border = c("', paste(col.border,collapse="\",\""), '"),\t # color of the borders
           col.bg = c("', paste(col.bg    ,collapse="\",\""), '"),\t # color of the backgrounds
          col.arr = c("', paste(col.arr   ,collapse="\",\""), '"),\t # color of the arrows
          lwd.arr = c(', paste( lwd.arr   ,collapse="," ),'),
         font.arr = c(', paste( font.arr  ,collapse="," ),'),
      col.txt.arr = c("', paste(col.txt.arr,collapse="\",\""), '"),\t # color of the text on arrows
       offset.arr =', offset.arr, ')\t\t # offset of parallel arrows \n',
 file=file )
}
invisible( list( Boxes = data.frame( pl.nam = pl.nam,
                                         xx = xx,
                                         yy = yy,
                                         wd = wd,
                                         ht = ht,
                                       font = font,
                                        lwd = lwd,
                                    col.txt = col.txt,
                                 col.border = col.border,
                                     col.bg = col.bg ),
                  Tmat = tm ) )
}
