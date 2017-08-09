######################################################################
### estimation method
LCa.fit <-
function( data, A, P, D, Y,
         model = "APa",    # or one of "ACa", "APaC", "APCa" or "APaCa" 
         a.ref,            # age reference for the interactions
        pi.ref = a.ref,    # age reference for the period interaction
        ci.ref = a.ref,    # age reference for the cohort interaction
         p.ref,            # period reference for the intercation
         c.ref,            # cohort reference for the interactions
          npar = c(a = 6,  # no. knots for main age-effect
                   p = 6,  # no. knots for peroid-effect
                   c = 6,  # no. knots for cohort-effect
                  pi = 6,  # no. knots for age in the period interaction
                  ci = 6), # no. knots for age in the cohort interaction
            VC = TRUE,     # numerical calculation of the Hessia?
         alpha = 0.05,     # 1 minus confidence level
           eps = 1e-6,     # convergence criterion
         maxit = 100,      # max. no iterations
         quiet = TRUE )    # cut the crap
{
# "model" must have values in c(APa/ACa/APaC/APCa/APaCa)?
if( !(model %in% c("APa","ACa","APaC","APCa","APaCa")) )
  stop( '"model" must be one of "APa","ACa","APaC","APCa","APaCa", but is', model,'\n' )

# Which main effects and interactions are in the model
 intP <- as.logical(length(grep("Pa",model)))
 intC <- as.logical(length(grep("Ca",model)))
mainP <- as.logical(length(grep("P" ,model))) # Also includes the age-period interaction
mainC <- as.logical(length(grep("C" ,model))) # Also includes the age-cohort product

# if a dataframe is supplied, fish out data and put in the function's environment
if( !missing(data) )
  {
  if (length(match(c("A", "P", "D", "Y"), names(data))) != 4)
  stop("Data frame ", deparse(substitute(data)),
       " has columns:\n", names(data),
       "\nmust have variables:\n", "A (age), P (period), D (cases) and Y (person-time)")
  data <- data[,c("A","P","D","Y")]
  data <- data[complete.cases(data),]
  A <- data$A
  P <- data$P
  D <- data$D
  Y <- data$Y
  } else { # if single vectors supplied, check they are all there
  nm <- c(missing(A),
          missing(P),
          missing(D),
          missing(Y))
  if (any(nm))
      stop("Variable", if (sum(nm) > 1)
          "s", paste(c(" A", " P", " D", " Y")[nm], collapse = ","),
          " missing from input")
  # and that they have the same length
  if( diff(range( lv <- c( length(A),
                           length(P),
                           length(D),
                           length(Y) ) )) != 0 )
      stop( "\nLengths of variables (", paste(paste(names(lv),
            lv, sep = ":"), collapse = ", "), ") are not the same." )
  } # end of data acquisition

# code-simplifier for knot calculation
eqqnt <- function(n) round( (1:n-0.5)/n, 2 )
# Define knots - we compute also the knots not needed
if( is.list(npar) ) {
  # Check if names is a named list
  if( is.null(names(npar)) ) stop( "If npar= is a list, it must be a *named* list.\n" )
    a.kn <- if( length(npar$a )>1 ) npar$a  else quantile( rep(  A,D), probs=eqqnt(npar$a ) )
    p.kn <- if( length(npar$p )>1 ) npar$p  else quantile( rep(P  ,D), probs=eqqnt(npar$p ) ) 
    c.kn <- if( length(npar$c )>1 ) npar$c  else quantile( rep(P-A,D), probs=eqqnt(npar$c ) )
   pi.kn <- if( length(npar$pi)>1 ) npar$pi else quantile( rep(  A,D), probs=eqqnt(npar$pi) )
   ci.kn <- if( length(npar$ci)>1 ) npar$ci else quantile( rep(  A,D), probs=eqqnt(npar$ci) )
    }
  else { # if npar is too short fill it up
  npar <- rep( npar, 5 )[1:5]
  # if not named, name it and notify
  if( is.null(names(npar)) ) names(npar) <- c("a","p","c","pi","ci")
    a.kn <- quantile( rep(  A,D), probs=eqqnt(npar["a"] ) )
    p.kn <- quantile( rep(P  ,D), probs=eqqnt(npar["p"] ) )
    c.kn <- quantile( rep(P-A,D), probs=eqqnt(npar["c"] ) )
   pi.kn <- quantile( rep(  A,D), probs=eqqnt(npar["pi"]) )
   ci.kn <- quantile( rep(  A,D), probs=eqqnt(npar["ci"]) )
  }
    
# Reference points
if( missing( p.ref) )  p.ref <- median( rep(P  ,D) )
if( missing( c.ref) )  c.ref <- median( rep(P-A,D) )
if( missing(pi.ref) ) pi.ref <- median( rep(  A,D) )    
if( missing(ci.ref) ) ci.ref <- median( rep(  A,D) )    

############################################################################
# Here starts the actual modelling
commence <- Sys.time()
    
# Matrices to extract the age-interaction terms at reference points
Mp <- Ns( rep(pi.ref,length(A)), knots=pi.kn, intercept=TRUE )
Mc <- Ns( rep(ci.ref,length(A)), knots=ci.kn, intercept=TRUE )

# Current age-effects (in the iteration these will be term predictions)
ba <- cbind( rep(1,length(A)), 1 )
# set to 0 if term is not in model at all
if( !mainP ) ba[,1] <- 0
if( !mainC ) ba[,2] <- 0

# Main effects model with (at least one) age-interactions
mat <- glm( D ~ -1 + Ns(   A, knots=a.kn, intercept=TRUE ) + 
                     Ns( P  , knots=p.kn, ref=p.ref):ba[,1] + 
                     Ns( P-A, knots=c.kn, ref=c.ref):ba[,2],
            offset = log(Y),
            family = poisson )
oldmb <- oldmat <- mat$deviance
    
# Terms prediction --- three terms here.
# No need to divide by the ba, it is eiter 1 or 0
pat <- predict( mat, type="terms" )

# iteration counter and continuation indicator
nit <- 0    
one.more <- TRUE
    
# For simple formatting of the iteration output
fC <- function(x,d) formatC(x,format="f",digits=d)
    
# now, iterate till convergence
while( one.more )
  {  
nit <- nit+1

# Terms with main effects should be either in interaction or offset,
# so one of these should always be 0
Pint <- Poff <- pat[,2]
Cint <- Coff <- pat[,3]
if( intP ) Poff <- Poff*0 else Pint <- Pint*0
if( intC ) Coff <- Coff*0 else Cint <- Cint*0
# Iteration of the age-interaction    
mb <- glm( D ~ -1 + Ns( A, knots=pi.kn, intercept=TRUE ):Pint +
                    Ns( A, knots=ci.kn, intercept=TRUE ):Cint,
           offset = pat[,1] + Poff + Coff + log(Y),
           family = poisson )

# Get the age-interaction terms only, and if one is not needed set to 0
ba <- predict( mb, type="terms" ) / cbind(Pint,Cint) /
      cbind( ci.lin( mb, subset="pi.kn", ctr.mat=Mp)[,1],
             ci.lin( mb, subset="ci.kn", ctr.mat=Mc)[,1] )
ba[is.na(ba)] <- 0

# If no interaction only main should be fitted; if no main effect, set to 0 
if( !intP ) ba[,1] <- rep(1,length(A)) * mainP
if( !intC ) ba[,2] <- rep(1,length(A)) * mainC
# apc model with assumed known interactions with age    
mat <- glm( D ~ -1 + Ns(   A, knots=a.kn, intercept=TRUE ) +
                     Ns( P  , knots=p.kn, ref=p.ref):ba[,1] + 
                     Ns( P-A, knots=c.kn, ref=c.ref):ba[,2],
            offset = log(Y),
            family = poisson )

# extract age and period terms      
pat <- predict( mat, type="terms" ) / cbind( 1, ba )
pat[is.na(pat)] <- 0

# convergence? Check bot that the two models give the same deviance
# and that the chnage in each is small
newmat <- mat$deviance
newmb  <-  mb$deviance
conv <- ( reldif <- max( (abs(newmat-newmb)/
                             (newmat+newmb)*2),
                         (oldmat-newmat)/newmat,
                         (oldmb -newmb )/newmb ) ) < eps
one.more <- ( !conv & ( nit < maxit ) )
oldmat <- newmat
oldmb  <- newmb
if( !quiet & nit==1 )
    cat( "    Deviances: model(AT) model(A) Rel. diff.\n" )
if( !quiet ) cat( "Iteration", formatC( nit, width=3, flag=" "), "",
                  fC(mat$deviance,3),
                  fC( mb$deviance,3),
                  fC( reldif, 7 ), "\n" )
  } # end of iteration loop

# Deviance and d.f - there is a "+1" because the intercept is in both models
# but not explicit, (both models fitted with "-1"), hence the df.null
# is the total no. observations 
dev <- mb$deviance
df  <- mat$df.null - ( mb$df.null- mb$df.res # no. parms in mb
                    + mat$df.null-mat$df.res # no. parms in mat
                    - 1 )                    # common intercept
if(  conv ) cat( "LCa.fit convergence in ", nit,
                " iterations, deviance:", dev, "on", df, "d.f.\n")
if( !conv ) cat( "LCa.fit *not* converged in ", nit,
                 " iterations:\ndeviance (AT):", mat$deviance,
                              ", deviance (B):" , mb$deviance, "\n",
                 if( VC ) "...no variance-covariance computed.\n" )
fin <- Sys.time()
if( !quiet ) cat("...using", round(difftime(fin,commence,units="secs"),1), "seconds.\n")
    
# unique values of A, P and C in the dataset for reporting effects
a.pt <- sort(unique(  A))
p.pt <- sort(unique(P  ))
c.pt <- sort(unique(P-A))

# extract effects from final models after convergence
ax <- ci.exp( mat, subset= "a.kn", ctr.mat=Ns(a.pt,knots= a.kn,intercept=TRUE ) )
kp <- ci.exp( mat, subset= "p.kn", ctr.mat=Ns(p.pt,knots= p.kn,ref=p.ref) )
kc <- ci.exp( mat, subset= "c.kn", ctr.mat=Ns(c.pt,knots= c.kn,ref=c.ref) )
pi <- ci.exp( mb , subset="pi.kn", ctr.mat=Ns(a.pt,knots=pi.kn,intercept=TRUE), Exp=FALSE )
ci <- ci.exp( mb , subset="ci.kn", ctr.mat=Ns(a.pt,knots=ci.kn,intercept=TRUE), Exp=FALSE )

# Label the estimated effects
rownames( ax ) <-    
rownames( pi ) <-    
rownames( ci ) <- a.pt
rownames( kp ) <- p.pt
rownames( kc ) <- c.pt

# do we bother about the correct variance-covariance?
if( VC & conv ) # ...certainly not without convergence
{
commence <- Sys.time()
if( !quiet ) cat("...computing Hessian by numerical differentiation...\n")

# the number of parameters for each of the 5 effects
na  <- length( grep(  "a.kn", names(coef(mat)) ) )
np  <- length( grep(  "p.kn", names(coef(mat)) ) )
nc  <- length( grep(  "c.kn", names(coef(mat)) ) )
npi <- length( grep( "pi.kn", names(coef(mb )) ) )
nci <- length( grep( "ci.kn", names(coef(mb )) ) )

# get only the parameters for effects that are non-zero (the others
# are in the models but they are 0)
ml.cf <- c( coef(mat)[c(rep(TRUE,na), 
                        rep(mainP,np),
                        rep(mainC,nc))],
             coef(mb)[c(rep(intP,npi),
                        rep(intC,nci))] )
# and some more snappy names for the parameters: first all names
all.nam <- c( paste("ax",1:na,sep=""),
              paste("kp",1:np,sep=""),
              paste("kc",1:nc,sep=""),
              paste("pi",1:npi,sep=""),
              paste("ci",1:nci,sep="") )
# ...then those actually present in the model
names( ml.cf ) <- all.nam[c(rep( TRUE,na),
                            rep(mainP,np),
                            rep(mainC,nc),
                            rep(intP,npi),
                            rep(intC,nci))]

# We need the variance-covariance of the estimates as the 2nd
# derivative of the log-likelihood, D*log(lambda) - lambda*Y,
# or for eta=log(lambda), D*eta - exp(eta)*Y,
# assuming the sequence of parameters is ax, kp, kc, pi, ci
# (first A, P, C from model mat, then Pa, Ca from model mb)
# Note that we cannot simplify this calculation because the model is
# non-linear in pi,kp resp. ci,kc

# Matrices to use in calculation of the terms of the model for each parms
MA  <- Ns(   A, knots= a.kn, intercept=TRUE )
Mp  <- Ns( P  , knots= p.kn, ref=p.ref )     
Mc  <- Ns( P-A, knots= c.kn, ref=c.ref )     
Mpi <- Ns(   A, knots=pi.kn, intercept=TRUE )
Mci <- Ns(   A, knots=ci.kn, intercept=TRUE )

# Computing the log-likelihood for any set of parameters
llik <-
function( parms )
{
              ax <- MA  %*% parms[   1:na]  ; nn <- na
if( mainP ) { kp <- Mp  %*% parms[nn+1:np]  ; nn <- nn+np } else kp = rep(0,length(ax))
if( mainC ) { kc <- Mc  %*% parms[nn+1:nc]  ; nn <- nn+nc } else kc = rep(0,length(ax))
if( intP  ) { pi <- Mpi %*% parms[nn+1:npi] ; nn <- nn+npi} else pi = rep(1,length(ax))
if( intC  ) { ci <- Mci %*% parms[nn+1:nci]               } else ci = rep(1,length(ax))
eta <- ax + pi*kp + ci*kc
sum( D*eta - exp(eta)*Y )
}

# Numerical calculation of the Hessian for llik
ivar <- -numDeriv::hessian( llik, ml.cf )

# Sometimes not quite positive definite, fix that after inverting the Hessian
vcov <- Matrix::nearPD( solve( ivar ) )
vcov <- as.matrix( vcov$mat )
fin <- Sys.time()
if( !quiet ) cat("...done - in", round(difftime(fin,commence,units="secs"),1), "seconds.\n")
    
# Since we now have the variances of the parameters for each of the
# effects we can compute corrected c.i.s for the effects
se.eff <-
function( sub, cM, Alpha=alpha )
{
wh <- grep( sub, names(ml.cf) )
res <- cbind(             cM %*% ml.cf[wh],
              sqrt( diag( cM %*% vcov[wh,wh] %*% t(cM) ) ) ) %*% ci.mat(alpha=Alpha)
colnames(res)[1] <- paste( "Joint", colnames(res)[1] )
res
}

# Append the corrected c.i.s to the effect objects
            ax <- cbind( ax, exp( se.eff( "ax", Ns(a.pt,knots= a.kn,intercept=TRUE) ) ) )
if( mainP ) kp <- cbind( kp, exp( se.eff( "kp", Ns(p.pt,knots= p.kn,ref=p.ref) ) ) )
if( mainC ) kc <- cbind( kc, exp( se.eff( "kc", Ns(c.pt,knots= c.kn,ref=c.ref) ) ) )
if( intP  ) pi <- cbind( pi,      se.eff( "pi", Ns(a.pt,knots=pi.kn,intercept=TRUE) ) )
if( intC  ) ci <- cbind( ci,      se.eff( "ci", Ns(a.pt,knots=ci.kn,intercept=TRUE) ) )
}

# Collect refs and knots in lists
klist <- list( a.kn=a.kn, pi.kn=pi.kn, p.kn=p.kn, ci.kn=ci.kn, c.kn=c.kn )
rlist <- list( pi.ref=pi.ref, p.ref=p.ref, ci.ref=ci.ref, c.ref=c.ref )
    
# Finally output object
res <- list( model = model,
                ax = ax, 
                pi = if( intP  ) pi else NULL, 
                kp = if( mainP ) kp else NULL, 
                ci = if( intC  ) ci else NULL, 
                kc = if( mainC ) kc else NULL, 
            mod.at = mat, 
            mod.b  = mb,
              coef = if( VC & conv ) ml.cf else NULL,
              vcov = if( VC & conv ) vcov  else NULL,
             knots = klist,
              refs = rlist,  
          deviance = dev,
       df.residual = df,
              iter = nit )
# Remove redundant stuff before returning
res <- res[!sapply(res,is.null)]
class( res ) <- "LCa"
invisible( res ) 
}

######################################################################
### utility to determine the types of terms in a model
model.terms <- 
function( x ) list( intP = as.logical(length(grep("Pa",x$model))),                
                   mainP = as.logical(length(grep("P" ,x$model))),                
                    intC = as.logical(length(grep("Ca",x$model))),                
                   mainC = as.logical(length(grep("C" ,x$model))) )

######################################################################
### print method
print.LCa <-
function( x, ... )
{
# terms in the model
mt <- model.terms( x )
    
# no. of parameters in each term    
npar <- sapply(x$knots,length)-1
npar[1] <- npar[1] + 1
npar <- npar[c(TRUE,mt$intP,mt$mainP,mt$intC,mt$mainC)]
npar <- paste( paste( npar, c(rep(", ",length(npar)-2)," and ",""), sep=""), collapse=rep("") )
cat( paste( x$model, ": Lee-Carter model with natural splines:\n",
     "  log(Rate) = ax(Age)",
     if( mt$mainP ) " + ",
     if(  mt$intP ) "pi(Age)",
     if( mt$mainP ) "kp(Per)",
     if( mt$mainC ) " + ",
     if(  mt$intC ) "ci(Age)",
     if( mt$mainC ) "kc(Coh)", "\nwith ",
     npar, " parameters respectively.\n",
     "Deviance: ", round(x$deviance,3), " on ", x$df, " d.f.\n",
     sep=""),
     if( !("vcov" %in% names(x)) ) "(only conditional c.i.s for effects)\n" )
}

######################################################################
### summary method
summary.LCa <-
function( object, show.est=FALSE, ... )
{
# terms in the model
mt <- model.terms( object )

print( object )

# the number of parameters for each of the 5 effects
na  <- length( grep(  "a.kn", names(coef(object$mod.at)) ) )
np  <- length( grep(  "p.kn", names(coef(object$mod.at)) ) )
nc  <- length( grep(  "c.kn", names(coef(object$mod.at)) ) )
npi <- length( grep( "pi.kn", names(coef(object$mod.b )) ) )
nci <- length( grep( "ci.kn", names(coef(object$mod.b )) ) )

# Show knots used    
cat( "\nKnots used:\n")
for( nm in names(object$knots[c(TRUE,mt$mainP,mt$intP,mt$mainC,mt$intC)]) )
   { cat( nm,"\n" ) ; print(object$knots[[nm]] ) }
cf <- rbind( ci.lin(object$mod.at)[c(rep(TRUE    ,na),   
                                     rep(mt$mainP,np),  
                                     rep(mt$mainC,nc)),1:2],
             ci.lin(object$mod.b )[c(rep(mt$intP ,npi),  
                                     rep(mt$intC ,nci)),1:2] )
if( "vcov" %in% names(object) )
  {  
  cf <- cbind( cf, sqrt( diag(object$vcov) ) )
  colnames( cf )[-1] <- c("cond.se","joint.se")
  }
if( show.est ) print( cf )
invisible( cf )
}

######################################################################
### plot method
plot.LCa <-
function( x, ... )
{
# terms in the model
mt <- model.terms( x )

# A small plot utility to exploit the structure of the effects    
plc <- 
function( x, ... ) matplot( as.numeric( rownames(x) ), x[,ncol(x)-2:0],
                            type="l", lty=1, lwd=c(3,1,1), ... )

plc( x$ax, col="black", xlab="Age", ylab="Age-specific rates", log="y" )
rug( x$knots$a.kn, lwd=2 )

if( mt$intP ){    
plc( x$pi, col="black", xlab="Age", ylab="Relative period log-effect multiplier" )
abline(h=1,v=x$refs$pi.ref)
rug( x$knots$pi, lwd=2 )
}

if( mt$mainP ){
plc( x$kp, col="black", log="y", xlab="Date of follow-up (period)", ylab="Period effect (RR)" )
abline(h=1,v=x$refs$p.ref)
rug( x$knots$kp, lwd=2 )
}

if( mt$intC ){    
plc( x$ci, col="black", xlab="Age", ylab="Relative cohort log-effect multiplier" )
abline(h=1,v=x$refs$ci.ref)
rug( x$knots$ci, lwd=2 )
}

if( mt$mainC ){    
plc( x$kc, col="black", log="y", xlab="Date of birth (cohort)", ylab="Cohort effect (RR)" )
abline(h=1,v=x$refs$c.ref)
rug( x$knots$kc, lwd=2 )
}

}

######################################################################
### predict method
predict.LCa <-
function( object,
         newdata,
           alpha = 0.05,
           level = 1-alpha,
             sim = ( "vcov" %in% names(object) ),
             ... )
{
# What main effects and interactions are in the model    
mt <- model.terms( object )    

# is person-years suppied, otherwise use units as in the model
if( "Y" %in% names(newdata) ) Y <-            newdata$Y else
                              Y <- rep(1,nrow(newdata))

# Matrices to extract effects at newdata rows
 Ma <- Ns(           newdata$A, knots = object$knots$a.kn, intercept = TRUE)
 Mp <- Ns( newdata$P          , knots = object$knots$p.kn, ref=object$refs$p.ref )
 Mc <- Ns( newdata$P-newdata$A, knots = object$knots$c.kn, ref=object$refs$c.ref )
Mpi <- Ns(           newdata$A, knots = object$knots$pi.kn, intercept = TRUE)
Mci <- Ns(           newdata$A, knots = object$knots$ci.kn, intercept = TRUE)

# Default terms values for models without interactions
kp <- kc <- rep( 0, nrow(newdata) )
pi <- ci <- rep( 1, nrow(newdata) )

# P, C and interaction term(s) if included in the model
if( mt$intP ) {
  pi <- ci.lin( object$mod.b , subset="pi.kn", ctr.mat=Mpi )[,1]
  kp <- ci.lin( object$mod.at, subset= "p.kn", ctr.mat=Mp  )[,1]
  }
if( mt$intC ) {
  ci <- ci.lin( object$mod.b , subset="ci.kn", ctr.mat=Mci )[,1]
  kc <- ci.lin( object$mod.at, subset= "c.kn", ctr.mat=Mc  )[,1]
  }

# First fitted values from mod.at
# Note that the model object mod.at always has the same number of
# parameters, for some of the models either period or cohort parameters
# are 0, hence not used.  
pr0 <- ci.exp( object$mod.at, alpha=alpha, ctr.mat=cbind(Ma,Mp*pi,Mc*ci) )

# Then fitted values from mod.b
# But mod.b has an offset beyond log(Y), namely all the APC terms
lp.b <-                ci.lin( object$mod.b , ctr.mat=cbind(Mpi*kp,Mci*kc) )[,1:2]
lp.b[,1] <- lp.b[,1] + ci.lin( object$mod.at, ctr.mat=cbind(Ma,Mp*(!mt$intP),Mc*(!mt$intC)) )[,1]
pr0 <- cbind( pr0, exp( lp.b %*% ci.mat(alpha=alpha) ) )

# label the estimates
colnames( pr0 )[c(1,4)] <- c("at|b Est.","b|at Est.")

# The doings above gives confidence intervals based on the conditional
# models, so if we want proper intervals we should simulate instead,
# using the posterior distribuion of all parameters, albeit under the
# slightly fishy assumption that the joint posterior is normal...
if( sim ) # also renders TRUE if sim is numerical (and not 0)
  {
if( is.logical(sim) & sim ) sim <- 1000
# Check that there is a vcov component of the model
if( !( "vcov" %in% names(object) ) )
    warning(
    "No variance-covariance in LCa object, only conditional c.i.s available.\n",
    "no simulation (prametric bootstrap) is done.\n" )  
else {   
# require( MASS )    
# using the parametric bootstrap based on the parameters and the
# (numerically computed) Hessian
eta <- NArray( list( pt = 1:nrow(pr0),
                     it = 1:sim ) )                      
parms <- MASS::mvrnorm( n = sim, 
                       mu = object$coef,
                    Sigma = object$vcov )
na  <- ncol( Ma  )
npi <- ncol( Mpi )
 np <- ncol( Mp  )
nci <- ncol( Mci )
 nc <- ncol( Mc  )
# Compute the linear predictor in each of the simulated samples
# period and cohort effects if not in the model
kp <- kc <- rep( 0, nrow(newdata) )
pi <- ci <- rep( 1, nrow(newdata) )
for( i in 1:sim ){
                 ax <- Ma  %*% parms[i,   1:na]  ; nn <- na
if( mt$mainP ) { kp <- Mp  %*% parms[i,nn+1:np]  ; nn <- nn+np }
if( mt$mainC ) { kc <- Mc  %*% parms[i,nn+1:nc]  ; nn <- nn+nc }
if( mt$intP  ) { pi <- Mpi %*% parms[i,nn+1:npi] ; nn <- nn+npi}
if( mt$intC  ) { ci <- Mci %*% parms[i,nn+1:nci]               }
eta[,i] <- ax + kp*pi + kc*ci
                  }
# predicted rates with bootstrap confidence limits              
pr.sim <- exp( t( apply( eta, 1, quantile,
                         probs=c(0.5,alpha/2,1-alpha/2), 
                         na.rm=TRUE ) ) )
colnames( pr.sim )[1] <- "Joint est."
return( pr.sim )
  }
}
else return( pr0 )    
}


