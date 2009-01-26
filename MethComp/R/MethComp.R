MethComp <-
function( data,
        random = c("mi","ir"),
     intercept = TRUE,
         slope = intercept,
      n.chains = 4,
        n.iter = 2000,
      n.burnin = n.iter/2,
        n.thin = ceiling((n.iter-n.burnin)/1000),
bugs.directory = getOption("bugs.directory"),
         debug = FALSE,
bugs.code.file = "model.txt",
       clearWD = TRUE,
        bugsWD = "bugsWD",
     code.only = FALSE,
      ini.mult = 2,
           org = FALSE,
       program = "BRugs",
           ... )
{
# Check the availability of required packages
Got.coda <- require( coda )
Got.R2WB <- require( R2WinBUGS )
Got.BRugs<- require( BRugs )
if( !Got.coda | ! Got.R2WB | !Got.BRugs )
  stop( "Using the MethComp function requires that\n",
        "the packages 'R2WinBUGS' as well as 'coda' are installed.\n",
        "In addition WinBUGS or the R-package Brugs is required too.\n",
        "All installed packages are shown if you type 'library()'." )

if(  is.null(bugs.directory) &&
    !is.null(bugs.dir <- getOption("R2WinBUGS.bugs.directory")) )
bugs.directory <- bugs.dir
   
# Check that a dataframe is supplied
if( !is.data.frame(data) | missing( data ) )
stop( "A dataframe should be supplied as the first argument." )

# Check that data has item, method and repl
rq.nam <- c("meth","item","repl","y")
if( sum( !is.na( wh <- match( rq.nam, names( data ) ) ) ) < 4 ) stop(
"\nThe supplied dataframe misses column(s) named ", rq.nam[is.na(wh)], ".\n" )

# Make program choice case-insensitive
program <- tolower( program )
program <- if( program %in% c("brugs","openbugs","ob") ) "openbugs"
    else { if( program %in% c(         "winbugs","wb") )  "winbugs"
    else { if( program %in% c(      "jags","jag","jg") )     "jags"
           else
           stop( "\n\nProgram '", program, "' not supported!" )
         } }

# Is the location of WinBUGS supplied if needed ?
if( !code.only & is.null( bugs.directory ) & program=="winbugs" ) stop(
"\nYou must supply the name of the folder where WinBUGS is installed,",
"\neither by using the parameter bugs.directory=...,",
"\n    or by setting options(bugs.directory=...).",
"\nThe latter will last you for the rest of your session.\n" )

if( !Got.BRugs & program=="openbugs" )
  stop( "Using the MethComp function with Brugs/openbugs option requires",
        "that the BRugs package is installed\n" )

# Logical variables determining the type of model
MI <- ( "MI" %in% toupper( random ) )
IR <- ( "IR" %in% toupper( random ) )

# Make the table of replicates for printing before method names are wiped
TT <- tab.repl( data )

# Get the names of methods that actually appear in data IN THAT ORDER.
# This quirk is necessary to avoid accidental reordering of method names
# if they are not alphabetically ordered
meth.names <- names( tt <- table( data$meth ) )[tt>0]

# Convert the variables to numerical 1,2,3,... for the sake of BUGS
# It MUST be factor() NOT as.factor()
data$meth <- as.integer( factor( data$meth ) )
data$item <- as.integer( factor( data$item ) )
data$repl <- as.integer( factor( data$repl ) )

# Quantities needed later
N  <- nrow( data )
Nm <- length( unique( data$meth ) )
Ni <- length( unique( data$item ) )
# The replicate numbers overall
Rn <- sort( unique( data$repl ) )
Nr <- max( Rn )
# Number of replicates per (item,method)
Nrep <- with( data, max(apply(table(item,meth,repl),1:2,function(x)sum(x>0))) )

# Print an explanatory text of what goes on:
cat( "\nComparison of", Nm, "methods, using", N, "measurements",
     "\non", Ni, "items, with up to", Nrep, "replicate measurements,",
     "\n(replicate values are in the set:", Rn, ")",
     "\n(",Nm,"*",Ni,"*",Nrep,"=",Ni*Nm*Nrep,"):",
     "\n\nNo. items with measurements on each method:\n" )
print( tab.repl( data ) )
cat( if( code.only ) "\nBugs program for a model with"
     else "\nSimulation run of a model with",
     if( !intercept & !slope ) "\n- no bias (intercept==0, slope==1)",
     if(  intercept & !slope ) "\n- fixed bias (slope==1)",
     if( !intercept &  slope ) "\n- proportional bias (intercept==0)",
     if( !MI & !IR ) "\n- no random interactions:",
     if(  MI & !IR ) "\n- method by item interaction:",
     if( !MI &  IR ) "\n- item by replicate interaction:",
     if(  MI &  IR ) "\n- method by item and item by replicate interaction:",
     if( code.only & !missing( bugs.code.file) )
         paste( "is written to the file",
                paste( getwd(), bugs.code.file, sep="/" ), "\n" )
     else if ( !code.only )
          paste( "\n- using", n.chains, "chains run for",
                                n.iter, "iterations \n  (of which",
                              n.burnin, "are burn-in),",
                 if( n.thin==1 ) "\n- monitoring all values of the chain:"
                 else paste( "\n- monitoring every", n.thin, "values of the chain:" ),
                 "\n- giving a posterior sample of",
                 round( n.chains*(n.iter-n.burnin)/n.thin ), "observations.\n\n" )
     )
# Make sure that it is printed before WinBUGS is fired up
flush.console()

# Make up name for a directory for all the junk files to be used by
# WinBUGS and read.bugs
shell( paste("mkdir",bugsWD) )

# Compute the range of the y's, and expand it to the range
# used for the "true" values for each item and for sd's
u.range <- range( data$y ) + c(-1,1) * 0.1 * diff( range( data$y ) )

# Write the BUGS gode to a file (or optionally the screen)
write.bugs.code( intercept=intercept, slope=slope, MI=MI, IR=IR,
                 N=N, Nm=Nm, Ni=Ni, Nr=Nr, u.range=u.range,
                 file = if( code.only )
                          {
                          if( missing(bugs.code.file) ) ""
                          else bugs.code.file
                          }
                        else paste(bugsWD,bugs.code.file,sep="/") )

# Generate the appropriate list of inits for the chains:
list.ini <- make.inits( data=data, Nm=Nm,
                        intercept=intercept, slope=slope, MI=MI, IR=IR,
                        n.chains=n.chains, ini.mult=ini.mult )

# If we want to execute the BUGS code --- well, then get on with it:
if( !code.only )
{
# Construct the data input data to WinBUGS
data.list <- c( list( N=N, Ni=Ni, Nm=Nm ),
                if( IR ) list( Nr=Nr  ),
                as.list( data[,c("meth","item",if( IR )"repl","y")] ) )

# Run bugs
res <- bugs( data = data.list,
            param = names( list.ini[[1]] ),
            inits = list.ini,
       model.file = bugs.code.file,
         n.chains = n.chains,
           n.iter = n.iter,
         n.burnin = n.burnin,
           n.thin = n.thin,
         bugs.dir = bugs.directory,
      working.dir = bugsWD,
            debug = debug,
          program = program,
          codaPkg = TRUE )

# and read the result into an mcmc.list object
# --- different approach for WinBUGS and OpenBUGE
if( program == "winbugs"  )
  res <- read.bugs( res, quiet=TRUE )
if( program == "openbugs" )
  res <- sims.array.2.mcmc.list( res$sims.array )

# Having read all the data, we don't need the junk directory any more
if( clearWD ) shell( paste("rmdir /S/Q", bugsWD ) )

# Now produce a mcmc object with the relevant parameters

# First add dummy colums for alpha and beta if they are not in the model.
# This facilitates all subsequent calculations
if( !intercept )
{
alphas <- rbind( rep( 0, Nm ) )
colnames( alphas ) <- paste( "alpha[", 1:Nm, "]", sep="" )
res <- addcols.mcmc( res, alphas )
}
if( !slope )
{
betas <- rbind( rep( 1, Nm ) )
colnames( betas ) <- paste( "beta[", 1:Nm, "]", sep="" )
res <- addcols.mcmc( res, betas )
}

# Construct a new mcmc object with the translation parameters and variance
# components as columns.
new.res <- trans.mcmc( res, MI, IR, Nm = Nm,
                            meth.names = meth.names,
                              n.chains = n.chains )

# Return the mcmc.list of relevant parameters
MCobj <- new.res

# Give class and attributes to the resulting object
class( MCobj ) <- c( "MethComp", class( MCobj ) )
attr( MCobj, "random" )  <- random
attr( MCobj, "methods" ) <- meth.names
attr( MCobj, "data" )    <- data
attr( MCobj, "mcmc.par" )<- list( n.chains = n.chains,
                                    n.iter = n.iter,
                                  n.burnin = n.burnin,
                                    n.thin = n.thin,
                                       dim = dim(as.matrix(MCobj)) )
if( org ) attr( MCobj, "orginal" ) <- res

invisible( MCobj )
}
# In case only the bugs code was wanted, we return the inits:
else {
   cat("\n\nThe generated inits for the BUGS code are:\n\n")
   print( list.ini )
   invisible( list.ini )
     }
}

################################################################################
### addcols.mcmc
################################################################################
addcols.mcmc <-
function( obj, cols )
{
if( !inherits( obj, "mcmc" ) &&
    !inherits( obj, "mcmc.list" ) ) stop( "obj must be a mcmc(.list) object" )
if( inherits( obj, "mcmc.list" ) )
  return( as.mcmc.list( lapply( obj, addcols.mcmc, cols ) ) )
else
{
if( length(dim(cols))==1 && length(cols)<1 )
  cols <- cols[rep(1:nrow(cols),nrow(obj))[1:nrow(obj)]]
if( length(dim(cols))>1 && nrow(cols)<nrow(obj) )
  cols <- cols[rep(1:nrow(cols),nrow(obj))[1:nrow(obj)],]
xx <- cbind( obj, cols )
attr( xx, "mcpar") <- attr( obj, "mcpar")
class( xx ) <- "mcmc"
return( xx )
}
}

################################################################################
### make.inits
################################################################################
make.inits <-
function( data, Nm, intercept, slope, MI, IR, n.chains, ini.mult=2 )
{
# function to clean up a list
rm.null <- function( lst ) lst[!sapply( lst, is.null )]
# n.chains must be at least 2:
if( n.chains < 2 ) stop( "n.chains must be at least 2, it is ", n.chains )
# Get variance component estimates and construct inits
vcm <- BA.est( data, linked=IR, MI=MI )
SDs <- vcm$Var.comp[,1]
nsd <- names( SDs )
if( MI ) sig.mi  <- SDs[grep("MxI",nsd)][1:(Nm-(Nm==2))] # Only one if Nm==2
if( IR ) sig.ir  <- SDs[grep("IxR",nsd)]
         sig.res <- SDs[grep("resid",nsd)]

# Produce a list of lists of length n.chains;
# the first one starting "correctly",
# and it has to be a list with one list as the first element
list.ini <- list( rm.null(
            list( "alpha" = if( intercept ) vcm$Bias,
                   "beta" = if( slope  ) rep(1,Nm),
               "sigma.mi" = if( MI ) sig.mi,
               "sigma.ir" = if( IR ) sig.ir,
              "sigma.res" = sig.res ) ) )
# and the subsequent with perturbed starting values for the variance components
for( j in 2:n.chains )
{
list.add <- rm.null(
            list( "alpha" = vcm$Bias,
                   "beta" = rep(1,Nm),
               "sigma.mi" = if( MI )
                            {
                            if( Nm > 2 )
                            sig.mi * ini.mult^sample(-1:1,Nm,replace=TRUE)
                            else
                            sig.mi * ini.mult^sample(-1:1, 1,replace=TRUE)
                            },
               "sigma.ir" = if( IR )
                            sig.ir * ini.mult^sample(-1:1, 1,replace=TRUE),
              "sigma.res" = sig.res* ini.mult^sample(-1:1,Nm,replace=TRUE) ) )
list.ini <- c( list.ini, list( list.add ) )
}
list.ini
}

################################################################################
### write.bugs.code
################################################################################
write.bugs.code <-
function( intercept, slope, MI, IR, # Defining structure of model
          u.range, # The range of the true values is relevant for the
                   # application of the initial values of the variance
                   # components and of initial values of(alpha,beta)=(0,1)
          N, Nm, Ni, Nr, # Defining structure of the data.
          file="" )
{
#  Write the BUGS code according to the arguments
cat( "model {
      for(j in 1:Ni)
        {
        u[j] ~ dunif(", paste( u.range, collapse="," ), ")
        }
      for (i in 1:N)
        {
        y[i] ~ dnorm( mu[i],tau.res[meth[i]] )
        mu[i] <- ",
if( intercept ) "alpha[meth[i]] +",
if( slope ) "
                  beta[meth[i]] *", "( u[item[i]]",
if( IR ) "+ e.ir[item[i],repl[i]]",
if( MI ) "+
                  e.mi[meth[i],item[i]]",
         ")",
        "
        }",
if( IR ) "
        for(r in 1:Nr)
           {
           for(i in 1:Ni)  { e.ir[i,r] ~ dnorm( 0, tau.ir ) }
           }",
if( MI ) paste( "
        for(m in 1:Nm)
           {
           for(i in 1:Ni)  { e.mi[m,i] ~ dnorm( 0,", if( Nm>2 ) "tau.mi[m]" else
                                                                "tau.mi", " ) }
           }" ),
if( IR ) paste("
                sigma.ir ~ dunif(0,", ceiling(10*u.range[2]), ") ;  tau.ir <- pow( sigma.ir,-2)" ),
if( MI & Nm == 2 )
      paste("
             sigma.mi ~ dunif(0,", ceiling(10*u.range[2]), ") ; tau.mi  <- pow(sigma.mi,-2)" ),
      "
       for( m in 1:Nm )
          { ",
if( intercept ) "
                alpha[m] ~ dnorm(0,0.0025)",
if( slope ) "
                 beta[m] ~ dunif(0,10)",
if( MI & Nm > 2 )
      paste("
             sigma.mi[m] ~ dunif(0,", ceiling(10*u.range[2]), ") ; tau.mi[m]  <- pow(sigma.mi[m],-2)" ),
       "
            sigma.res[m] ~ dunif(0,", paste( ceiling(10*u.range[2]) ), ") ; tau.res[m] <- pow(sigma.res[m],-2)
          }
}
       ",
       file = file )
}

################################################################################
### abconv
################################################################################
abconv <-
function( a1, b1=1:4, a2=NULL, b2=NULL,
          col.names=c("alpha.2.1","beta.2.1","id.2.1") )
{
if( ( inherits( a1, "data.frame" ) |
      inherits( a1, "matrix" ) )
    & length( b1 )==4 )
  {
  cols <- a1
  wh <- b1
  a1 <- cols[,wh[1]]
  a2 <- cols[,wh[2]]
  b1 <- cols[,wh[3]]
  b2 <- cols[,wh[4]]
  }
a2.1 <- a2 - a1 * b2 / b1
b2.1 <- b2 / b1
id2.1 <- a2.1 / ( 1-b2.1 )
dfr <- data.frame( a2.1, b2.1, id2.1 )
names( dfr ) <- col.names
dfr
}

################################################################################
### conv.par
################################################################################
conv.par <-
function( sim.mat, i, j )
{
# Function to produce posteriors of the conversion parameters from method i to j.

# Find the columns of sim.mat with the relevant alphas and betas
# Note the "\\" necessary to escape the special meaning of "[" in grep.
wh <-  c( grep( paste("alpha\\[",i,"]",sep=""), colnames(sim.mat) ),
          grep( paste("alpha\\[",j,"]",sep=""), colnames(sim.mat) ),
          grep( paste( "beta\\[",i,"]",sep=""), colnames(sim.mat) ),
          grep( paste( "beta\\[",j,"]",sep=""), colnames(sim.mat) ) )

# First compute the alpha, beta and intersection with the identity
# The result is a matrix object with 3 columns
ab.conv <- abconv( sim.mat, wh,
                   col.names = paste( c("alpha","beta","id"), j, i, sep="." ) )
                     
# The variance for predicting method j from i:
# First the residual variances
var.conv <- sim.mat[,paste("sigma.res[",j,"]", sep="" )]^2 +
            ab.conv[,paste("beta",j,i,sep=".")]^2 *
            sim.mat[,paste("sigma.res[",i,"]", sep="" )]^2

# And if there is a method by item interaction this one too
# Recall that the sigma.mi random effect is specified as multiplied by beta
MI <- any( as.logical( grep( "sigma.mi\\[", colnames(sim.mat) ) ) )
if( MI ) var.conv <- var.conv +
                     sim.mat[,paste("beta[", j, "]", sep="" )]^2 *
                     sim.mat[,paste("sigma.mi[", j, "]", sep="" )]^2 +
                     ab.conv[,paste("beta",j,i,sep=".")]^2 *
                     sim.mat[,paste("beta[", i, "]", sep="" )]^2 *
                     sim.mat[,paste("sigma.mi[", i, "]", sep="" )]^2
sd.conv <- sqrt( var.conv )
ab.conv <- cbind( ab.conv, sd.conv )
names( ab.conv ) = paste( c("alpha","beta","id.int","sd.pred"),
                          "[", j, "].[", i, "]", sep="" )
ab.conv <- as.matrix( ab.conv )
if( i == j ) ab.conv[,4,drop=FALSE]
        else ab.conv
}

################################################################################
### sims.array.2.mcmc.list
################################################################################
sims.array.2.mcmc.list <-
function( aa )
{
zz <- list(list())
for( i in 1:(dim(aa)[2]) )
   {
   tmp <- coda:::mcmc( aa[,i,] )
   zz <- c( zz, list(tmp) )
   }
return( coda:::mcmc.list( zz[-1] ) )
}

################################################################################
### mat2.mcmc
################################################################################
mat2.mcmc.list <-
function( mm, n.chains )
{
zz <- list(list())
n.sims <- dim(mm)[1]/n.chains
if( floor(n.sims)!=n.sims)
  stop( "Matrix supplied does not have nrows a multiple of n.chains" )
for( i in 1:n.chains )
{
tmp <- mcmc( mm[(i-1)*n.sims+(1:n.sims),] )
zz <- c( zz, list(tmp) )
}
return( mcmc.list( zz[-1] ) )
}

################################################################################
### trans.mcmc
################################################################################
trans.mcmc <-
function( res, MI, IR, names=TRUE, Nm, meth.names, n.chains )
{
# First convert to a matrix
sim.mat <- as.matrix( res )

# Get the residual variances 'as is'
new.res <- sim.mat[,grep("sigma.res",colnames(sim.mat))]

# Parameters for each method
for( i in Nm:1 )
{
# Conversion from i to j
for( j in Nm:1 )
{
new.res <- cbind( conv.par(sim.mat,i,j), new.res )
}
# The remaining variance components:
# Put the sds of the MI and IR effects on the right scales:
if( MI )
  {
  # First if Nm==2 there is only 1 sigma.mi which must be multiplied with
  # each of the two estimated betas
  if( Nm==2 )
    {
    new.res <- cbind( new.res,
                      sim.mat[,"sigma.mi"]*
                      sim.mat[,paste("beta[",i,"]",sep="")] )
    }
  else
    {
    new.res <- cbind( new.res,
                      sim.mat[,paste("sigma.mi[",i,"]",sep="")]*
                      sim.mat[,paste(    "beta[",i,"]",sep="")] )
    }
  colnames( new.res )[ncol(new.res)] <- paste("sigma.mi[",i,"]",sep="")
  }
if( IR )
  {
  new.res <- cbind( new.res,
                    sim.mat[,"sigma.ir"]*
                    sim.mat[,paste("beta[",i,"]",sep="")] )
  colnames( new.res )[ncol(new.res)] <- paste("sigma.ir[",i,"]",sep="")
  }
}

# Total variance for each method (as SD, of course)
for(i in 1:Nm )
  {
  var.tot <- new.res[,paste("sigma.res[",i,"]",sep="")]^2
  if( IR ) var.tot <- var.tot + new.res[,paste("sigma.ir[",i,"]",sep="")]^2
  if( MI ) var.tot <- var.tot + new.res[,paste("sigma.mi[",i,"]",sep="")]^2
  new.res <- cbind( new.res, sqrt( var.tot ) )
  colnames( new.res )[ncol(new.res)] <- paste("sigma.tot[",i,"]",sep="")
  }

# Put the method names into the posterior results matrix colnames
if( names )
{
zz <<- colnames( new.res )
for( i in 1:Nm )
   zz <- gsub( paste("\\[",i,"]",sep=""),
               paste("[",meth.names[i],"]",sep=""), zz )
zz <- gsub( "].\\[", ".", zz )
zz -> colnames( new.res )
}

return( mat2.mcmc.list( new.res, n.chains ) )
}
