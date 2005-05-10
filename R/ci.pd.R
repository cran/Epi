ci.pd <-
  function( aa, bb=NULL, cc=NULL, dd=NULL, alpha=0.05, print = TRUE )
  {
# Computes the c.i. for the probability difference as
# method 10 from Newcombe, Stat.Med. 17, (1998), pp.873 ff.

# Allow various forms of matrix input    
  if( is.vector( aa ) & length( aa ) > 1 ) print <- FALSE 
  if( length( dim( aa ) ) == 2 )
    {
    bb <- aa[1,2]
    cc <- aa[2,1]
    dd <- aa[2,2]
    aa <- aa[1,1]
    }
  if( length( dim( aa ) ) == 3 )
    {
    bb <- aa[1,2,]
    cc <- aa[2,1,]
    dd <- aa[2,2,]
    aa <- aa[1,1,]
    print <- FALSE
    }
  if( length( dim( aa ) ) > 3 )
    stop( "Maximal array dimension (3) exceeded!" )

# Function to give roots in a 2nd degree polynomial
# (Polyroot does not work on vectors of coefficients)
  pol2 <-
  function( Aye, Bee, Sea )
  {
  Dee <- Bee^2 - 4 * Aye * Sea
  lo <- ifelse( Dee >= 0, ( -Bee - sqrt( Dee ) ) / ( 2 * Aye ), NA )
  hi <- ifelse( Dee >= 0, ( -Bee + sqrt( Dee ) ) / ( 2 * Aye ), NA )
  cbind( lo, hi )
  }
  if( print ) twoby2( rbind( c(aa,cc), c(bb,dd) ) ) 
# Put the data in the right form  
  x1 <- aa
  n1 <- aa+cc
  x2 <- bb
  n2 <- bb+dd
  zz <- qchisq( 1-alpha, 1 )
  A1 <-        1 + zz / n1
  B1 <- -2*x1/n1 - zz / n1
  C1 <-          ( x1 / n1 )^2
  r1 <- pol2( A1, B1, C1 )
  A2 <-        1 + zz / n2
  B2 <- -2*x2/n2 - zz / n2
  C2 <-          ( x2 / n2 )^2
  r2 <- pol2( A2, B2, C2 )
  pd <- x1/n1 - x2/n2
  dlt <- sqrt( ( x1/n1 - r1[1] )^2 +
               ( x2/n2 - r2[2] )^2 )
  eps <- sqrt( ( x1/n1 - r1[2] )^2 +
               ( x2/n2 - r2[1] )^2 )             
  res <- cbind( pd, pd-dlt, pd+eps )
  colnames( res ) <- c("prob.diff","lo","hi")
  rownames( res ) <- paste( aa, "/(", aa, "+", cc, ") - ",
                            bb, "/(", bb, "+", dd, ") = ", sep="" )
  invisible( res )
  }
