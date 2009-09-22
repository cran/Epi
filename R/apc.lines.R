apc.lines <-
function( A,
          P,
          C,
      scale = c("log","ln","rates","inc","RR"),
  frame.par = NULL,
      drift = 0,
         c0 = median( C[,1] ),
         a0 = median( A[,1] ),
         p0 = c0 + a0,
         ci = rep( FALSE, 3 ),
        lwd = c(3,1,1),
        lty = 1,
        col = "black",
       type = "l",
      knots = FALSE,
        ...
        )
{
# What scale are we using:
log.scale <- ( match.arg( scale ) %in% c("log","ln") )
# Are confidence intervals requested ?
# Expand ci to three components if only one is given
if( is.logical( ci ) ) ci <- rep( ci, 3 )[1:3]
# Allow ci to be any character string or vector with "a", "p" or "c" in it
if( is.character( ci ) )
  {
  ci <- paste( ci, collapse="" )
  ci <- unlist( strsplit( toupper( ci ), split="" ) )
  ci <- !is.na( match( c("A","P","C"), ci ) )
  }
# Check the format of the various inputs are correct
if( inherits( A, "apc" ) )
  {
 obj <- A
  p0 <- ifelse( is.na( A$Ref[1] ), median( A$Per[,1] ), A$Ref[1] )
  c0 <- ifelse( is.na( A$Ref[2] ), median( A$Coh[,1] ), A$Ref[2] )
  a0 <- p0 - c0
   C <- A$Coh
   P <- A$Per
   A <- A$Age
  log.scale <- FALSE
  }
if( ( ci[1] & dim( A )[2] < 4 ) | dim( A )[2] < 2 )
  stop( "A, ", deparse( substitute( A ) ), ") must be a ", 2+2*ci[1]," column matrix" )
if( ( ci[2] & dim( P )[2] < 4 ) | dim( P )[2] < 2 )
  stop( "P, ", deparse( substitute( P ) ), ") must be a ", 2+2*ci[2]," column matrix" )
if( ( ci[3] & dim( C )[2] < 4 ) | dim( C )[2] < 2 )
  stop( "C, ", deparse( substitute( C ) ), ") must be a ", 2+2*ci[3]," column matrix" )
if( p0 != c0+a0 )
  stop( "p0=", p0, " differs from c0 + a0 = ", c0, " + ", a0, " = ", c0+a0 )
# Expand single graphics parameters to length 3
  col  <- rep(  col, 3 )
  lty  <- rep(  lty, 3 )
  type <- rep( type, 3 )
# Transform to log-scale if input is as rates and RRs
  if( !log.scale )
    {
    if( missing( drift ) ) drift <- 1
      drift <- log( drift )
    A[,2:4] <- log( A[,2:4] )
    P[,2:4] <- log( P[,2:4] )
    C[,2:4] <- log( C[,2:4] )
    }
  A[,2:4] <- exp( A[,2:4] - drift * ( A[,1] - a0 ) )
  P[,2:4] <- exp( P[,2:4] + drift * ( P[,1] - p0 ) )
  C[,2:4] <- exp( C[,2:4] - drift * ( C[,1] - c0 ) )
  # If no frame parameters are given
  if( is.null( frame.par ) )
  frame.par <- c( min( C[,1] ) + max( A[,1] ),
                  ifelse( scale=="log",
                          exp( mean(      A[,2] ) ),
                          exp( mean( log( A[,2] ) ) ) ) )
  # Now we can plot the lines
  matlines( A[,1], A[,ifelse( ci[1], -1, 2)],
            col=col[1], lwd=lwd, lty=lty[1], type=type[1], ... )
  matlines( P[,1] - frame.par[1], P[,ifelse( ci[2], -1, 2)] * frame.par[2],
            col=col[2], lwd=lwd, lty=lty[2], type=type[2], ... )
  matlines( C[,1] - frame.par[1], C[,ifelse( ci[3], -1, 2)] * frame.par[2],
            col=col[3], lwd=lwd, lty=lty[3], type=type[3], ... )
  points( obj$Ref - frame.par[1], frame.par[c(2,2)], pch=16, col="white" )
  points( obj$Ref - frame.par[1], frame.par[c(2,2)], pch=1, lwd=2, col=col[2:3] )
  if( knots & inherits( obj, "apc" ) ) 
  { 
  rug( obj$Knots$Age, side=1, col=col[1] )
  rug( obj$Knots$Per - frame.par[1], side=1, col=col[2] )
  rug( obj$Knots$Coh - frame.par[1], side=3, col=col[3] )
  }             
  box()
  invisible( list( A=A, P=P, C=C ) )
}
