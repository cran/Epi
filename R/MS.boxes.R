tbox <-
function( txt, x, y, w, h,
          font=2, txt.col="black",
          lwd=2, border="black", col="transparent" )
{
rect( x-w/2, y-h/2, x+w/2, y+h/2, lwd=lwd, border=border, col=col )
text( x, y, txt, font=font, col=txt.col )
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
function( x1, y1, x2, y2, fr=0.8,
          angle=17, lwd=2, length=par("pin")[1]/30, ... )
{
if( fr > 1 ) fr <- fr/100
for( a in 1:angle )
arrows( x1 + (x2-x1)*(1-fr)/2,
        y1 + (y2-y1)*(1-fr)/2,
        x2 - (x2-x1)*(1-fr)/2,
        y2 - (y2-y1)*(1-fr)/2, angle=a, lwd=lwd, ... )
}

boxarr <-
function( b1, b2, ... )
{
x1 <- b1[1]
y1 <- b1[2]
w1 <- b1[3]
h1 <- b1[4]
x2 <- b2[1]
y2 <- b2[2]
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
invisible( list( x=(bx1+bx2)/2, y=(by1+by2)/2 ) )
}

boxes <- function (obj, ...) UseMethod("boxes")

boxes.Lexis <-
function( obj, file, mult=1.5, cex=1.5 )
{
if( inherits(obj,"Lexis") ) tm <- tmat( obj )
else if( is.matrix(obj) & diff(dim(obj))==0 ) tm <- obj
else stop( "First argument must be a Lexis object or a square matrix.\n" )

# State names
st.nam <- colnames( tm )
if( is.null(st.nam) ) st.nam <- paste(1:ncol(tm))
n.st <- length( st.nam )

# Here is the plot
par( mar=c(0,0,0,0), cex=cex )
plot( NA,
      bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )
xx <- yy <- wd <- ht <- numeric(n.st)
b <- list()
for( i in 1:n.st )
 {
 txt <- st.nam[i]
 cat( "\nClick for level ", txt )
 flush.console()
 pt <- locator(1)
 xx[i] <- pt$x
 yy[i] <- pt$y
 ht[i] <- strheight( txt ) * mult * 1.5
 wd[i] <- strwidth( txt ) * mult
 b[[i]] <- tbox( txt, xx[i], yy[i], wd[i], ht[i] )
 }
cat( "\n" )

for( i in 1:n.st ) for( j in 1:n.st )
  if( !is.na(tm[i,j]) )
    boxarr( b[[i]], b[[j]] )
# If we want the code
if( !missing(file) )
{
xx <- round(xx,1)
yy <- round(yy,1)
wd <- round(wd,1)
ht <- round(ht,1)
obj.name <- substitute( deparse( obj ) )
cat( '
st.nam <- colnames( tm )
if( is.null(st.nam) ) st.nam <- paste(1:ncol(tm))
n.st <- length( st.nam )
b <- list(NULL)
par( mar=c(0,0,0,0), cex=1.5 )
plot( NA,
      bty="n",
      xlim=0:1*100, ylim=0:1*100, xaxt="n", yaxt="n", xlab="", ylab="" )',
      file = file )
for( i in 1:n.st )
cat( 'b[[',paste(i),']] <-\ntbox("',st.nam[i],'",',
     paste(c(xx[i],yy[i],wd[i],ht[i]),collapse=","),",",
'\n         font=2, txt.col="black",',
'\n         lwd=2, border="black", col="transparent" )\n',
     sep="", file=file, append=TRUE )
for( i in 1:n.st ) for( j in 1:n.st )
  if( !is.na(tm[i,j]) )
cat( 'boxarr( b[[',paste(i),']], b[[',paste(j),']], fr=0.8,
             angle=15, lwd=2, length=par("pin")[1]/30 )\n',
      sep="", file=file, append=T )
}
}