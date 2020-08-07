# The coef() methods in nlme and lme4 do something different,
# other objects do not even have coef or vcov methods defined,
# so we make a workaround by specifying our own generic methods:
COEF          <- function( x, ... ) UseMethod("COEF")
COEF.default  <- function( x, ... ) coef( x, ... )
VCOV          <- function( x, ... ) UseMethod("VCOV")
VCOV.default  <- function( x, ... ) vcov( x, complete=TRUE, ... )

# Then we can get from these methods what we want from lme, mer etc.
COEF.lme      <- function( x, ... ) nlme::fixed.effects( x )
COEF.mer      <- function( x, ... ) lme4::fixef( x )
COEF.lmerMod  <- function( x, ... ) lme4::fixef( x )
# The vcov returns a matrix with the wrong class so we strip that:
VCOV.lme      <- function( x, ... ) as.matrix(vcov( x ))
VCOV.mer      <- function( x, ... ) as.matrix(vcov( x ))
VCOV.lmerMod  <- function( x, ... ) as.matrix(vcov( x ))

# For the rest of the non-conforming classes we then just need the methods not defined
# VCOV.coxph    <- function( object, ... ) survival:A::vcov.coxph( object, complete=FALSE, ... )
COEF.crr      <- function( object, ... ) object$coef
VCOV.crr      <- function( object, ... ) object$var
COEF.MIresult <- function( object, ... ) object$coefficients
VCOV.MIresult <- function( object, ... ) object$variance
COEF.mipo     <- function( object, ... ) object$qbar
VCOV.mipo     <- function( object, ... ) object$t
COEF.polr     <- function( object, ... ) summary(object)$coefficients
VCOV.gnlm     <- function( object, ... ) object$cov
VCOV.rq       <- function( object, ... ) summary(object, cov=TRUE)$cov

df2ctr <-
function( obj, nd )
    {
if( !( inherits(obj,"lm") | inherits(obj,"coxph") ) )
  stop("data frame facility not inplemented for ", class(obj), " objects" )
# Factors in the prediction frame must have more than one level, which
# they typically do not have in the specification, so we find the
# factors in the prediction frame and expand levels to the complete
# set of levels which should secure the working of model.matrix()
dcl <- attr( obj$terms, "dataClasses" )
whf <- ( dcl == "factor" )
if( any(whf) ) for( fn in names(dcl)[which(whf)] )
                nd[,fn] <- factor( nd[,fn], levels=obj$xlevels[[fn]] )
# The contrast matrix from the model - differs a bit between (g)lm, gam and coxph
# this is needed to keep NA rows from the data frame supplied
org.op <- options( na.action='na.pass' )
on.exit( options( org.op ) )
if( inherits(obj,"coxph") ) MM <- model.matrix(         obj     ,    data=nd )
if( inherits(obj,"gam"  ) ) MM <- model.matrix(         obj     , newdata=nd )
    else
if( inherits(obj,"lm"   ) ) MM <- model.matrix( formula(obj)[-2],    data=nd )
return( MM )
    }

ci.dfr <-
function( obj, ndx, ndr,
         xvars = NULL,
          vcov = FALSE,
         alpha = 0.05,
           Exp = FALSE,
        sample = FALSE )
{
if( nrow(ndr)==1 ) ndr <- ndr[rep(1,nrow(ndx)),,drop=FALSE]
if( (    ( nrow(ndx) !=  nrow(ndr)) ) |
    ( any(names(ndx) != names(ndr)) ) )
    stop("\nThe two prediction frames must have same dimensions and column names:",
         "but dimensions are: (",
         paste( dim(ndx),collapse=","), ") and (",
         paste( dim(ndr),collapse=","), ")\n", 
         "and column names are:\n",
         "exp: ", paste( names(ndx), collapse=", " ), "\n",
         "ref: ", paste( names(ndr), collapse=", " ), "\n")
# Now supply and fix those variables that are needed in order to get model.matrix working:
# Supplied variable names:
 cols <- names( ndx )
# Factors in model; which are supplied; derive names of omitted factors (ofacs)
 facs <- names( obj$xlevels )
ofacs <- setdiff( facs, cols )
# omitted *variables* must be supplied
ovars <- setdiff( xvars, facs )
# Construct the extra columns with bogus data (their contribution will be null)
xcols <- ndx[,NULL]
if( length(ofacs) > 0 ) for( fn in ofacs ) xcols[,fn] <- obj$xlevels[[fn]][1]
if( length(ovars) > 0 ) for( vn in ovars ) xcols[,vn] <- 1
if( dim(xcols)[2]>0 )
  {
  ndx <- cbind( ndx, xcols )
  ndr <- cbind( ndr, xcols )
  }
ci.lin( obj,
    ctr.mat = df2ctr( obj, ndx ) - df2ctr( obj, ndr ),
       vcov = vcov,
      alpha = alpha,
        Exp = Exp,
     sample = sample )
}

ci.dfr2 <-
function( obj, nd1, nd2, nd3, nd4,
         xvars = NULL,
          vcov = FALSE,
         alpha = 0.05,
           Exp = FALSE,
        sample = FALSE )
{
if( nrow(nd2)==1 ) nd2 <- nd2[rep(1,nrow(nd1)),,drop=FALSE]
if( nrow(nd3)==1 ) nd3 <- nd3[rep(1,nrow(nd1)),,drop=FALSE]
if( nrow(nd4)==1 ) nd4 <- nd4[rep(1,nrow(nd1)),,drop=FALSE]
if( ( nrow(nd2) !=  nrow(nd1)) |
    ( nrow(nd3) !=  nrow(nd1)) |
    ( nrow(nd4) !=  nrow(nd1)) |
 any(names(nd2) != names(nd1)) |
 any(names(nd3) != names(nd1)) |
 any(names(nd4) != names(nd1)) )
    stop("\nThe prediction frames must have same dimensions and column names:",
         "but dimensions are: (",
         paste( dim(nd1),collapse=","), ") and (",
         paste( dim(nd2),collapse=","), ") and (",
         paste( dim(nd3),collapse=","), ") and (",
         paste( dim(nd4),collapse=","), ")\n", 
         "and column names are:\n",
         "1: ", paste( names(nd1), collapse=", " ), "\n",
         "2: ", paste( names(nd2), collapse=", " ), "\n",
         "3: ", paste( names(nd3), collapse=", " ), "\n",
         "4: ", paste( names(nd4), collapse=", " ), "\n")
# Now supply and fix those variables that are needed in order to get model.matrix working:
# Supplied variable names:
 cols <- names( nd1 )
# Factors in model; which are supplied; derive names of omitted factors (ofacs)
 facs <- names( obj$xlevels )
ofacs <- setdiff( facs, cols )
# omitted *variables* must be supplied
ovars <- setdiff( xvars, facs )
# Construct the extra columns with bogus data (their contribution will be null)
xcols <- nd1[,NULL]
if( length(ofacs) > 0 ) for( fn in ofacs ) xcols[,fn] <- obj$xlevels[[fn]][1]
if( length(ovars) > 0 ) for( vn in ovars ) xcols[,vn] <- 1
if( dim(xcols)[2]>0 )
  {
  nd1 <- cbind( nd1, xcols )
  nd2 <- cbind( nd2, xcols ) 
  nd3 <- cbind( nd3, xcols )
  nd4 <- cbind( nd4, xcols ) 
  }
ci.lin( obj,
    ctr.mat = df2ctr( obj, nd1 )
             -df2ctr( obj, nd2 )
             -df2ctr( obj, nd3 )
             +df2ctr( obj, nd4 ),
       vcov = vcov,
      alpha = alpha,
        Exp = Exp,
     sample = sample )
}

ci.lin <-
function( obj,
      ctr.mat = NULL,
       subset = NULL,
       subint = NULL,
        xvars = NULL,
        diffs = FALSE,
         fnam = !diffs,
         vcov = FALSE,
        alpha = 0.05,
           df = Inf,
          Exp = FALSE,
       sample = FALSE )
{
# If ctr.mat is a data frame, call df2ctr 
if( inherits( ctr.mat, "data.frame" ) ) ctr.mat <- df2ctr( obj, ctr.mat )

# If ctr.mat is a list of two dataframes then call ci.dfr
if( inherits( ctr.mat, "list" ) )
  {
  if( length(ctr.mat)==2 )
    {  
  if( !inherits( ctr.mat[[1]], "data.frame" ) |
      !inherits( ctr.mat[[2]], "data.frame" ) )
      stop("If ctr.mat is a list it must be a list of data frames")
  return( ci.dfr( obj, ctr.mat[[1]], ctr.mat[[2]],
                xvars = xvars,
                 vcov = vcov,
                alpha = alpha,
                  Exp = Exp,
               sample = sample ) )
    }
# If ctr.mat is a list of two dataframes then call ci.dfr2 for 2nd
# order differences
  if( length(ctr.mat)==4 )
    {  
  if( !inherits( ctr.mat[[1]], "data.frame" ) |
      !inherits( ctr.mat[[2]], "data.frame" ) |
      !inherits( ctr.mat[[3]], "data.frame" ) |
      !inherits( ctr.mat[[4]], "data.frame" ) )
      stop("If ctr.mat is a list it must be a list of data frames")
  return( ci.dfr2( obj, ctr.mat[[1]], ctr.mat[[2]], ctr.mat[[3]], ctr.mat[[4]],
                 xvars = xvars,
                  vcov = vcov,
                 alpha = alpha,
                   Exp = Exp,
                sample = sample ) )
    }
  }

# First extract all the coefficients and the variance-covariance matrix
cf  <- COEF( obj )
vcv <- VCOV( obj )
# Alised parameters are set to 0
if( any( is.na(  cf ) ) )  cf[is.na(cf )] <- 0
if( any( is.na( vcv ) ) ) vcv[is.na(vcv)] <- 0

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
}
# end of the function all.dif for all differences 

# Were all differences requested?
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
  # class( res ) <- c("ci.lin","matrix")
    }
# Return the requested structure
if( sample ) invisible( res )
else
if( vcov ) invisible( list( coef=ct[,1], vcov=vc ) ) else res
}

# print.ci.lin <-
# function( x, ..., digits=3 )
# {
# print( round( unclass(x), digits ) )
# }

# Handy wrapper
ci.exp <-
function( ..., Exp=TRUE, pval=FALSE )
{
res <- if( Exp ) 
         {
ci.lin( ..., Exp=TRUE  )[,if(pval) c(5:7,4)   else 5:7     ,drop=FALSE]
         } else {
ci.lin( ..., Exp=FALSE )[,if(pval) c(1,5,6,4) else c(1,5,6),drop=FALSE]
                }
# class( res ) <- c( "ci.lin", "matrix" )
res
}

# Wrapper for predict.glm to give estimates and confidence intervals
ci.pred <-
function( obj, newdata,
         Exp = NULL,
       alpha = 0.05 )
{
if( !inherits( obj, "lm" ) ) stop("Not usable for non-(g)lm objects")
if( !inherits( obj, "glm" ) & inherits( obj, "lm" ) )
    return( predict.lm( obj, newdata=newdata, interval="confidence" ) )
# get the prediction and se on the link scale
else {
zz <- predict( obj, newdata=newdata, se.fit=TRUE, type="link" )
# compute ci on link scale
zz <- cbind( zz$fit, zz$se.fit ) %*% ci.mat( alpha=alpha )
# transform as requested
if( missing(Exp) ) {   return( obj$family$linkinv(zz) )
} else {  if(  Exp ) { return(                exp(zz) ) 
   } else if( !Exp )   return(                    zz  )
       }
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
if( is.data.frame(r1) ) r1 <- as.matrix( r1 )
if( is.data.frame(r2) ) r2 <- as.matrix( r2 )
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
