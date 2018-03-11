# The coef() methods in nlme and lme4 do something different,
# other objects do not even hav coef of vcov methods defined,
# so we make a workaround by specifying our own generic methods:
COEF          <- function( x, ... ) UseMethod("COEF")
COEF.default  <- function( x, ... ) coef( x, ... )
VCOV          <- function( x, ... ) UseMethod("VCOV")
VCOV.default  <- function( x, ... ) vcov( x, complete=FALSE, ... )

# Then we can get from these methods what we want from lme, mer etc.
COEF.lme      <- function( x, ... ) nlme::fixed.effects( x )
COEF.mer      <- function( x, ... ) lme4::fixef( x )
COEF.lmerMod  <- function( x, ... ) lme4::fixef( x )
# The vcov returns a matrix with the wrong class so we strip that:
VCOV.lme      <- function( x, ... ) as.matrix(vcov( x ))
VCOV.mer      <- function( x, ... ) as.matrix(vcov( x ))
VCOV.lmerMod  <- function( x, ... ) as.matrix(vcov( x ))

# For the rest of the non-conforming classes we then just need the methods not defined
COEF.crr      <- function( object, ... ) object$coef
VCOV.crr      <- function( object, ... ) object$var
COEF.MIresult <- function( object, ... ) object$coefficients
VCOV.MIresult <- function( object, ... ) object$variance
COEF.mipo     <- function( object, ... ) object$qbar
VCOV.mipo     <- function( object, ... ) object$t
COEF.polr     <- function( object, ... ) summary(object)$coefficients
VCOV.gnlm     <- function( object, ... ) object$cov
VCOV.rq       <- function( object, ... ) summary(object, cov=TRUE)$cov

# Computing the contrasts between two prediction data frames
# --- a utility called from ci.lin 
ci.dfr <-
function( obj, ndx, ndr, Exp )
    {
if( !is.data.frame(ndx) ) stop("2nd argument in list must be a data frame")
if( !is.data.frame(ndr) ) stop("3rd argument in list must be a data frame")
if( nrow(ndr)==1 ) ndr <- ndr[rep(1,nrow(ndx)),,drop=FALSE]
if( (    ( nrow(ndx) !=  nrow(ndr)) ) |
    ( any(names(ndx) != names(ndr)) ) )
    stop("\nThe two prediction frames must have same dimensions and column names:",
         "but dimensions are: ", dim(ndx), " and ", dim(ndr), "\n", 
         "and column names are:\n",
         "exp: ", names(ndx), "\n",
         "ref: ", names(ndr), "\n")
# Now fix those variable that are needed in order to get model.matrix working
# Supplied variable names
 cols <- names( ndx )
# Factors in model; which are supplied; derive names of omitted
 facs <- names(obj$xlevels)
ofacs <- match( facs, cols )
ofacs <- facs[is.na(ofacs)]
# Variables in model; which are supplied; derive names of omitted
 vars <- setdiff( all.vars(obj$formula)[-1], facs )
ovars <- match( vars, cols )
ovars <- vars[is.na(ovars)]
# Construct the extra columns
xcols <- ndx[,NULL]
if( length(ofacs) > 0 ) for( fn in ofacs ) xcols <- cbind( xcols, obj$xlevels[[fn]][1] )
if( length(ovars) > 0 ) for( vn in ovars ) xcols <- cbind( xcols, 1 )
if( dim(xcols)[2]>0 )
  {
  names( xcols ) <- c(ofacs,ovars)
  ndx <- cbind( ndx, xcols )
  ndr <- cbind( ndr, xcols )
  }
# Factors must have more than one level, which they typically will not
# have in the specification of the prediction frames, so we find the
# factors in the model and expand levels to the complete set of levels
dcl <- attr(obj$terms,"dataClasses")
whf <- ( dcl == "factor" )
if( any(whf) )
  for( fn in names(dcl)[which(whf)] )
     {
     ndx[,fn] <- factor( ndx[,fn], levels=obj$xlevels[[fn]] )
     ndr[,fn] <- factor( ndr[,fn], levels=obj$xlevels[[fn]] )
     }
# Then we can set up the contrast matrix and call ci.lin
CM <- model.matrix( formula(obj)[-2], data=ndx ) -
      model.matrix( formula(obj)[-2], data=ndr )
ci.lin( obj, ctr.mat=CM, Exp=Exp )
    }

ci.lin <-
function( obj,
      ctr.mat = NULL,
       subset = NULL,
       subint = NULL,
        diffs = FALSE,
         fnam = !diffs,
         vcov = FALSE,
        alpha = 0.05,
           df = Inf,
          Exp = FALSE,
       sample = FALSE )
{
# If ctr.mat is a list of two dataframes call ci.dfr
if( is.list(ctr.mat) )
  {
  if( !is.data.frame(ctr.mat[[1]]) |
      !is.data.frame(ctr.mat[[2]]) )
      stop("If ctr.mat is a list it must be a list of two data frames")
  return( ci.dfr( obj, ctr.mat[[1]], ctr.mat[[2]], Exp=Exp ) )
  }
  
# First extract all the coefficients and the variance-covariance matrix
cf  <- COEF( obj )
vcv <- VCOV( obj )
# Workaround to expand the vcov matrix with 0s so that it matches
# the coefficients vector in case of (extrinsic) aliasing.
if( any( is.na( cf ) ) )
  {
  if( inherits( obj, c("coxph") ) )
    { # aliased parameters are only NAs in coef, but omitted from vcov
    wh <- !is.na(cf)
    cf <- cf[wh]
    vcv <- vcv[wh,wh]
    }
  else
    {  
  if( inherits( obj, c("clogistic") ) )
    {
    cf[is.na(cf)] <- 0
    }
  else
    {
    vM <- matrix( 0, length( cf ), length( cf ) )
    dimnames( vM ) <- list( names( cf ), names( cf ) )
    vM[!is.na(cf),!is.na(cf)] <- vcv
    # vM <- vcv
    cf[is.na(cf)] <- 0
    vcv <- vM
    }
    }
  }

# Function for computing a contrast matrix for all possible
# differences between a set of parameters.
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
# end of the function all.dif for all differences 
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
         vcv <- vcov( obj, complete=FALSE )[wh,wh]
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
  if( is.character( subint ) ) {
    sb <- 1:length(cf)
    for( i in 1:length(subint) ) sb <- intersect( sb, grep(subint[i],names(cf)) )
    subset <- sb # unique( sb )
    }
  if( is.null( subset ) & is.null( subint ) ) subset <- 1:length( cf )
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
  if( sample )
    {
    # mvrnorm() returns a vector if sample=1, otherwise a sample by
    # length(cf) matrix - hence the rbind so we always get a row 
    # matrix and res then becomes an nrow(ctr.mat) by sample matrix
    res <- ctr.mat %*% t( rbind(mvrnorm( sample, cf, vcv )) )   
    }  
  else
    {
    ct <- ctr.mat %*% cf
    vc <- ctr.mat %*% vcv %*% t( ctr.mat )
    se <- sqrt( diag( vc ) )
    ci <- cbind( ct, se ) %*% ci.mat( alpha=alpha, df=df )
    t0 <- cbind( se, ct/se, 2 * ( pnorm( -abs( ct / se ) ) ) )
    colnames(t0) <- c("StdErr", "z", "P")
    res <- cbind(ci, t0)[, c(1, 4:6, 2:3), drop=FALSE]
    if( Exp )
      {
      res <- cbind(      res[,1:4     ,drop=FALSE],
                    exp( res[,c(1,5,6),drop=FALSE] ) )
      colnames( res )[5] <- "exp(Est.)"
      }
    }
# Return the requested structure
if( sample ) invisible( res ) else
if( vcov ) invisible( list( coef=ct[,1], vcov=vc ) ) else res
}

# Handy wrapper
ci.exp <-
function( ..., Exp=TRUE, pval=FALSE )
{
if( Exp )
ci.lin( ..., Exp=TRUE  )[,if(pval) c(5:7,4)   else 5:7     ,drop=FALSE]
else
ci.lin( ..., Exp=FALSE )[,if(pval) c(1,5,6,4) else c(1,5,6),drop=FALSE]
}

# Wrapper for predict.glm to give estimates and confidnece intervals
ci.pred <-
function( obj, newdata,
         Exp = NULL,
       alpha = 0.05 )
{
if( !inherits( obj, "glm" ) ) stop("Not usable for non-glm objects")
# get the prediction and se on the link scale
zz <- predict( obj, newdata=newdata, se.fit=TRUE, type="link" )
# compute ci on link scale
zz <- cbind( zz$fit, zz$se.fit ) %*% ci.mat( alpha=alpha )
# transform as requested
if( missing(Exp) ) {   return( obj$family$linkinv(zz) )
} else {  if(  Exp ) { return(                exp(zz) ) 
   } else if( !Exp )   return(                    zz  )
       }  
}

# Function to calculate RR with CIs from independent rates with CIs;
# r1 and r2 are assumed to be vectors or 2 or 3-column matrices with
# rate, lower and upper confidence limits repectively.
ci.ratio <-
function( r1, r2,
         se1 = NULL, # standard error of rt1
         se2 = NULL, # standard error of rt2
      log.tr = !is.null(se1) & !is.null(se2), # is this log-rates?
       alpha = 0.05,
        pval = FALSE )
{
if( is.matrix(r1) & !is.null(se1) ) warning("r1 is matrix, se1 is ignored")
if( is.matrix(r2) & !is.null(se2) ) warning("r2 is matrix, se2 is ignored")

# if supplied as 1-column matrix change to vector
if( is.matrix(r1) ) if( ncol(r1)==1 ) r1 <- as.vector( r1 )
if( is.matrix(r2) ) if( ncol(r2)==1 ) r2 <- as.vector( r2 )

# move to log scale
if( !log.tr ) {
  r1 <- log( r1 ) 
  r2 <- log( r2 ) 
  }

# how wide are the condidence intervals    
if( is.matrix(r1) ) if( ncol(r1)>1 ) rg1 <- t( apply(r1,1,range) )
if( is.matrix(r2) ) if( ncol(r2)>1 ) rg2 <- t( apply(r2,1,range) )

# get the estimates on the log-scale
R1 <- if( is.matrix(r1) ) apply( rg1, 1, mean ) else r1 
R2 <- if( is.matrix(r2) ) apply( rg2, 1, mean ) else r2 
if( is.null(se1) ) se1 <- apply( rg1, 1, diff ) / (2*qnorm(1-alpha/2))
if( is.null(se2) ) se2 <- apply( rg2, 1, diff ) / (2*qnorm(1-alpha/2))

# compute the RR and the c.i. and optionally the p-value
 lrr <- R1 - R2
slrr <- sqrt( se1^2 + se2^2 )
  rr <- cbind(lrr,slrr) %*% ci.mat(alpha=alpha)
if( !log.tr ) rr <- exp( rr )
if( pval ) return( cbind( rr, 1-pchisq( (lrr/slrr)^2, 1 ) ) )
      else return(        rr                                )    
}
