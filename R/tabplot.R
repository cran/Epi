tabplot <-
function( M,
        col,
     border = "black",
        lwd = 2,
    collabs = TRUE,
    rowlabs = NULL,
      equal = FALSE,
        las = 1,
       main = NULL,
   cex.main = 1.0,
      vaxis = FALSE
         )
{
 # Number of rows and columns
 #
 nr <- nrow( M )
 nc <- ncol( M )

 # Coloring defaults to graytones
 #
 if(       missing( col ) ) col <- gray( nr:1 / nr ) else
 if(  is.function ( col ) ) col <- col( nr )         else
 if(  is.character( col ) ) col <- col[1:nr]         else
 if( !is.character( col ) ) stop(
     "col= argument wrong mode (should be function or character" )
 
 # Compute the relative size of the column totals
 #
 mg <- apply( M, 2, sum )
 bw <- cumsum( mg ) / sum ( mg )
 # But don't use it if equal is requested
 if( equal ) bw <- 1:nc / nc

 # Plot everything in a unit square
 #
 plot( 0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i" )

 # Calculate column and total percentages
 #
 pctM <- sweep( M, 2, mg, "/" )
 tpcM <- M / sum( M ) * 100

 # Matrices holding the corners of the rectangles
 #
 xl <- yb <- xr <- yt <- matrix( NA, nr, nc )
 for( i in 1:nc )
   {
      xl[,i] <- c(0,bw[-nc])[i]
      yb[,i] <- c(0,cumsum( pctM[,i] )[-nr])
      xr[,i] <- bw[i]
      yt[,i] <- cumsum( pctM[,i] )
    }
 
 # Then plot the rectangles with a border to make it nice
 #
 rect( xl, yb, xr, yt, col=rep(col,nc), border=border, lwd=lwd )
 box( lwd=lwd, col=border )

 # If column labels are needed put them in the right places
 #
 if( collabs ) axis( side=1, at = bw - diff( c( 0, bw ) ) / 2,
                     labels=colnames( M ), tick=FALSE, las=las )

 # Put row labels correctly at none, one or both of the sides
 #
 llab <- FALSE
 if( !is.null( rowlabs ) )
   {
 if( "L" %in% unlist( strsplit( toupper( rowlabs ), "" ) ) ){
   ccol <- cumsum( pctM[,1] )
   axis( side = 2, at = ccol - diff( c( 0, ccol ) ) / 2,
         labels=rownames( M ), tick=FALSE, las=las )
   llab <- TRUE
   }
 if( "R" %in% unlist( strsplit( toupper( rowlabs ), "" ) ) ){
   ccol <- cumsum( pctM[,nc] )
   axis( side = 4, at = ccol - diff( c( 0, ccol ) ) / 2,
         labels=rownames( M ), tick=FALSE, las=las )
   }
   }

 # If an axis is requested, put it on
 if( is.character( vaxis ) ){
 if( "L" %in% unlist( strsplit( toupper( vaxis ), "" ) ) ){
   axis( side = 2 )
   }
 if( "R" %in% unlist( strsplit( toupper( vaxis ), "" ) ) ){
   axis( side = 4 )
   }
   } else
 if( is.logical( vaxis ) & vaxis ){
   axis( side=ifelse( llab, 4, 2 ) )
   } else
 if( !is.logical( vaxis ) ) stop( "'vaxis' must be logical or  charcater." )
 
 
 # Title if requested
 #
 if( !is.null( main ) ) mtext( side=3, main, line=1, cex=cex.main )

 # Return the corners of the rectangles as a list
 #
 invisible( list( xl=xl, yb=yb, xr=xr, yt=yt, pct=tpcM ) )
}
