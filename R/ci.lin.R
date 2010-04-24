ci.lin <-
function( obj,
      ctr.mat = NULL,
       subset = NULL,
        diffs = FALSE,
         fnam = !diffs,
         vcov = FALSE,
        alpha = 0.05,
           df = Inf,
          Exp = FALSE ) 
{
# First extract all the coefficients and the variance-covariance matrix
#
if( any( inherits( obj, c("coxph","glm","gls","lm","nls","survreg","clogistic") ) ) ) {
       cf <- coef( obj )
      vcv <- vcov( obj )
} else if( inherits( obj, c("lme") ) ) {
       cf <- fixed.effects( obj )
      vcv <- vcov( obj )
} else if( inherits( obj, c("mer") ) ) {
       cf <- fixef( obj )
      vcv <- as.matrix( vcov( obj ) )
} else if (inherits(obj, "MIresult")) {
       cf <- obj$coefficients
      vcv <- obj$variance
} else if( inherits( obj, "polr" ) ) {
       cf <- summary( obj )$coefficients
      vcv <- vcov( obj )
} else if( inherits( obj, "gnlm" ) ) {
       cf <- coef( obj )
      vcv <- obj$cov
} else stop( "\"", deparse( substitute( obj ) ), "\" is of class \"",
              class( obj ), "\" which is not supported." )

# Workaround to expand the vcov matrix with 0s so that it matches
# the coefficients vector in case of (extrinsic) aliasing.
if( any( is.na( cf ) ) )
  {
if( inherits( obj, c("coxph") ) )
  {
  wh <- !is.na(cf)
  cf <- cf[wh]
  vcv <- vcv[wh,wh]
  }
else
if( inherits( obj, c("clogistic") ) )
  {
  cf[is.na(cf)] <- 0
  }
else
  {
vM <- matrix( 0, length( cf ), length( cf ) )
dimnames( vM ) <- list( names( cf ), names( cf ) )
vM[!is.na(cf),!is.na(cf)] <- vcv
cf[is.na(cf)] <- 0
vcv <- vM
  }
  }

# Function for computing a contrast matrix for all possible
# differences between a set of parameters.
#
all.dif <-
function( cf, pre=FALSE )
{
nn <- length( cf )
nr <- nn * ( nn - 1 ) / 2
nam <- names( cf )

# Work out the indexes of parameter pairs to compare
#
xx <- numeric( 0 )
for( j in 2:nn ) xx <- c(xx, j:nn )
ctr <- cbind( rep( 1:(nn-1), (nn-1):1 ), xx )

# Now for the annotation:
# Find out how large a proportion of rownames are identical
i <- 1
while( all( substr( nam, 1, i ) == substr( nam[1], 1, i ) ) ) i <- i+1 

# If a factor name is given, then use this, otherwise the identical part
# of the parameter names
if( is.character( pre ) ) {
prefix <- pre
  pre <- TRUE
} else {
prefix <- substr( nam[1], 1, i-1 )
}
rn <- paste( if( pre ) prefix else "",
             substring( nam[ctr[,1]], i ), "vs.",
             substring( nam[ctr[,2]], i ) )

# Finally, construct the contrast matrix and attach the rownames
cm <- matrix( 0, nr, nn )
cm[cbind(1:nr,ctr[,1])] <- 1 
cm[cbind(1:nr,ctr[,2])] <- -1 
rownames( cm ) <- rn
cm
}

# Were all differences requested?
#
if( diffs )
  {
  if( is.character( subset ) )
    {
    if ( inherits( obj, "lm" ) &
         length( grep( subset, names( obj$xlevels ) ) )>0 )
       { # The case of factor level differences we find the relevant
         # subset of parameters by reconstructing names of parameters
         wf <- grep( subset, af <- names( obj$xlevels ) )
         # All factor levels
         fn <- obj$xlevels[[af[wf]]]
         # Reconstruct names of relevant parameter names
         pnam <- paste( af[wf], fn, sep="" )
         # Find them in the parameter vector
         wh <- match( pnam, names( coef( obj ) ) )
         # Get the relevant subset, and stick in 0s for NAs
         cf <- coef( obj )[wh]
         cf[is.na( cf )] <- 0
         vcv <- vcov( obj )[wh,wh]
         vcv[is.na( vcv )] <- 0
         names( cf ) <- rownames( vcv ) <- colnames( vcv ) <-
             paste( subset, ": ", fn, sep="" )
       } else {
         subset <- grep( subset, names( cf ) )
             cf <-  cf[subset]
            vcv <- vcv[subset,subset] 
       }
    } else {
     cf <-  cf[subset]
    vcv <- vcv[subset,subset] 
   }
   ctr.mat <- all.dif( cf, pre=fnam )
  }

if( !diffs )
  {
  if( is.character( subset ) ) {
    sb <- numeric(0)
    for( i in 1:length( subset ) ) sb <- c(sb,grep( subset[i], names( cf )  ))
    subset <- sb # unique( sb )
    }
  if( is.null( subset ) ) subset <- 1:length( cf )
  # Exclude units where aliasing has produced NAs.
  # Not needed after replacement with 0s
  # subset <- subset[!is.na( cf[subset] )]
   cf <-  cf[subset]
  vcv <- vcv[subset,subset]
  if( is.null( ctr.mat ) )
    {
    ctr.mat <- diag( length( cf ) )
    rownames( ctr.mat ) <- names( cf )
    }
  if( dim( ctr.mat )[2] != length(cf) )
      stop( paste("\n Dimension of ", deparse(substitute(ctr.mat)), 
            ": ", paste(dim(ctr.mat), collapse = "x"),
            ", not compatible with no of parameters in ", 
            deparse(substitute(obj)), ": ", length(cf), sep = ""))
  }
# Finally, here is the actual computation
    ct <- ctr.mat %*% cf
    vc <- ctr.mat %*% vcv %*% t( ctr.mat )
    se <- sqrt( diag( vc ) )
    ci <- cbind( ct, se ) %*% ci.mat( alpha=alpha, df=df )
    t0 <- cbind( se, ct/se, 2 * ( 1 - pnorm( abs( ct / se ) ) ) )
    colnames(t0) <- c("StdErr", "z", "P")
    res <- cbind(ci, t0)[, c(1, 4:6, 2:3), drop=FALSE]
    if( Exp ) {
      res <- cbind(      res[,1:4     ,drop=FALSE],
                    exp( res[,c(1,5,6),drop=FALSE] ) )
      colnames( res )[5] <- "exp(Est.)"
    }
# Return the requested structure
if( vcov ) invisible( list( est=ct, vcov=vc ) ) else res
}
