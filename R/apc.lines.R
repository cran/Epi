apc.lines <-
function( A,
          P,
          C,
      scale = "log",
  frame.par = c( min( C[,1] ) + max( A[,1] ),
                 ifelse( scale=="log",
                         exp( mean(      A[,2] ) ),           
                         exp( mean( log( A[,2] ) ) ) ) ),
      drift = 0,
         c0 = median( C[,1] ),
         a0 = median( A[,1] ),
         p0 = c0 + a0,
         ci = rep( FALSE, 3 ),
        lwd = c(3,1,1),
        lty = rep( 1, 3),
        col = rep( "black", 3 )
        )
{
  # Check the format of the various inputs are correct
  if( ( ci[1] & dim( A )[2] < 4 ) | dim( A )[2] < 2 )
    stop( "A, ", deparse( substitute( A ) ), ") must be a ", 2+2*ci[1]," column matrix" )
  if( ( ci[2] & dim( P )[2] < 4 ) | dim( P )[2] < 2 )
    stop( "P, ", deparse( substitute( P ) ), ") must be a ", 2+2*ci[2]," column matrix" )
  if( ( ci[3] & dim( C )[2] < 4 ) | dim( C )[2] < 2 )
    stop( "C, ", deparse( substitute( C ) ), ") must be a ", 2+2*ci[3]," column matrix" )
  if( p0 != c0+a0 )
    stop( "p0=", p0, " differs from c0 + a0 = ", c0, " + ", a0, " = ", c0+a0 )
  # Allow only a single color to be given
  col <- rep( col, 3 )
  # What scale are we using:
  if( ! scale %in% c("log","ln","rates","inc","RR") )
    stop( 'scale must be one of "log", "ln", "rates", "inc", "RR"' )
  log.scale <- ( toupper( scale ) %in% c("LOG","LN") )
  # Are confidence intervals requested ?
  if( is.logical( ci ) ) ci <- rep( ci, 3 )
  # Just allow ci be any character string or vector with
  # "a","p" or "c" in it
  if( is.character( ci ) )
    {
    ci <- paste( ci, collapse="" )
    ci <- unlist( strsplit( toupper( ci ), split="" ) )
    ci <- !is.na( match( c("A","P","C"), ci ) )
    }
  # Transform to log-scale if input is as rates and RRs
  if( !log.scale )
    {
    drift <- log( drift )
    A[,2:4] <- log( A[,2:4] )
    P[,2:4] <- log( P[,2:4] )
    C[,2:4] <- log( C[,2:4] )
    }
  A[,2:4] <- exp( A[,2:4] - drift * ( A[,1] - a0 ) )
  P[,2:4] <- exp( P[,2:4] + drift * ( P[,1] - p0 ) )
  C[,2:4] <- exp( C[,2:4] - drift * ( C[,1] - c0 ) )
  matlines( A[,1], A[,2 + ci[1] * 0:2],
            col=col[1], lwd=lwd, lty=lty[1] )
  matlines( P[,1] - frame.par[1], P[,2 + ci[2] * 0:2] * frame.par[2],
            col=col[2], lwd=lwd, lty=lty[2] )
  matlines( C[,1] - frame.par[1], C[,2 + ci[3] * 0:2] * frame.par[2],
            col=col[3], lwd=lwd, lty=lty[3] )
  invisible( list( A=A, P=P, C=C ) )
}
