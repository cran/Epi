plotMeth <-
function( data,
         which = NULL,
        col.LA = "blue",
      cex.name = 2,
     var.range,
    diff.range,
     var.names = FALSE,
           ... )
{
if( "meth" %in% names( data ) ) data <- to.wide( data )
if( missing( which ) )
    which <- levels( attr( data, "reshapeWide" )$times )
if( is.null( which ) ) stop( "Please indicate variables to plot" )

# Should we plot the variable names
if( is.logical( var.names ) )
  {
  plot.names <- var.names
  if( is.character( which ) ) var.names <- which
  else  var.names <- names( data )[which]
  }
else plot.names <- TRUE

# Functions to plot the upper and lower panels
pnl <-
function( x, y, ... )
{
abline( 0, 1, col="black" )
points( x, y, pch=16, ... )
}

pnu <-
function( x, y, ... )
{
sdd <-   sd( (y-x), na.rm=TRUE )
mnd <- mean( (y-x), na.rm=TRUE )
abline( h=mnd+(-1:1)*2.00*sdd, col=col.LA )
abline( h=0,col="black" )
points( (x+y)/2, y-x, pch=16, ... )
}

pldat <- data[,which]
nvar <- ncol( pldat )
if( missing(var.range) )
  rg <- range( pldat, na.rm=TRUE )
  else rg <- var.range
if( missing(diff.range) )
  dif.rg <- c(-1,1)*diff(rg)/2
  else if( length( diff.range )==1 ) dif.rg <- c(-1,1)*abs(diff.range)
       else dif.rg <- diff.range

# Make a grid of plots preserving the old plotting paramters
oma <- par("mar")
oldpar <- par( mfrow=c(nvar,nvar), mar=c(0,0,0,0), oma=c(4,4,4,4), las=1 )
on.exit( par ( oldpar ) )
for( ir in 1:nvar ) for( ic in 1:nvar )
   {
   if( ir == ic )
     {
     plot( 0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE )
     text( 0.5, 0.5, names( pldat )[ir], font=2, cex=cex.name )
     box()
     }
   if( ir < ic )
     {
     plot( 0:1, 0:1, xlim=rg, ylim=dif.rg,
           type="n", xlab="", ylab="", axes=FALSE )
     pnu( pldat[,ir], pldat[,ic], ... )
     if( plot.names )
       {
       text( rg[1], dif.rg[2], paste(var.names[ic],"-",var.names[ir]), adj=c(0,1) )
       text( rg[2], dif.rg[1], paste("(",var.names[ic],"+",var.names[ir],")/2"), adj=c(1,0) )
       }
     if( ir == nvar ) axis( side=1 )
     if( ic == 1 )    axis( side=2 )
     if( ir == 1 )    axis( side=3 )
     if( ic == nvar ) axis( side=4 )
     box()
     }
   if( ir > ic )
     {
     plot( 0:1, 0:1, xlim=rg, ylim=rg,
           type="n", xlab="", ylab="", axes=FALSE )
     pnl( pldat[,ic], pldat[,ir], ... )
     if( plot.names )
       {
       text( rg[1], rg[2], var.names[ir], adj=c(0,1) )
       text( rg[2], rg[1], var.names[ic], adj=c(1,0) )
       }
     if( ir == nvar ) axis( side=1 )
     if( ic == 1 )    axis( side=2 )
     if( ir == 1 )    axis( side=3 )
     if( ic == nvar ) axis( side=4 )
     box()
     }
   }
}
