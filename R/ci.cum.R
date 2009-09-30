ci.cum <-
function( obj,
      ctr.mat = NULL,
       subset = NULL,
         intl = 1,
        alpha = 0.05,
          Exp = TRUE )
{
# First extract all the coefficients and the variance-covariance matrix
#
if( any( inherits( obj, c("coxph","glm","gls","lm","nls","survreg") ) ) ) {
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

# Check if the intervals matches ctr.mat
if( length( intl ) == 1 ) intl <- rep( intl, nrow( ctr.mat ) )
if( length( intl ) != nrow( ctr.mat ) ) stop( "intl must match ctr.mat" )

# Workaround to expand the vcov matrix with 0s so that it matches
# the coefficients vector in case of (extrinsic) aliasing.
if( any( is.na( cf ) ) )
  {
vM <- matrix( 0, length( cf ), length( cf ) )
dimnames( vM ) <- list( names( cf ), names( cf ) )
vM[!is.na(cf),!is.na(cf)] <- vcv
cf[is.na(cf)] <- 0
vcv <- vM
   }

# If subset is not given, make it the entire set
#
# if( is.null( subset ) ) subset <- 1:length( cf )

# Useful function for constructing a matrix linking estimate, s.e. to
# a confidence interval
ci.mat <-
function( alpha = 0.05 )
{
ciM <- rbind( c(1,1,1), qnorm(1-alpha/2)*c(0,-1,1) )
colnames( ciM ) <- c("Estimate",
   paste( formatC( 100*   alpha/2 , format="f", digits=1 ), "%", sep="" ),
   paste( formatC( 100*(1-alpha/2), format="f", digits=1 ), "%", sep="" ) )
ciM
}

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

# Finally, here is the actual computation of the estimates
    ct <- ctr.mat %*% cf
    vc <- ctr.mat %*% vcv %*% t( ctr.mat )
# If Exp was requested, we take the exponential of the estimates
# before we cumulate the sum
if( Exp )
{
    ct <- exp( ct )
    vc <- ( ct[,1] %*% t(ct[,1]) ) * vc
}
# Here is the cumulation matrix
    cum.mat <- 1 - upper.tri( diag(ct[,1]) )
    # Multiply columns of the matrix with interval lengths
    cum.mat <- t( intl * t( cum.mat ) )
# This is then multiplied to the coefficients
    ct <- cum.mat %*% ct
    vc <- cum.mat %*% vc %*% t( cum.mat )
    se <- sqrt( diag( vc ) )
    cbind( cbind( ct, se ) %*% ci.mat( alpha=alpha ), Std.err.=se )
}
