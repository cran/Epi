summary.MethComp <-
function( object, alpha=0.05,
          ... )
{
Got.coda <- require( coda )
if( !Got.coda )
  stop( "Using the summary.MethComp function requires that\n",
        "the package 'coda' is installed.\n",
        "All installed packages are shown if you type 'library()'." )

MI <- ( "mi" %in% attr( object, "random" ) )
mnam <- attr( object, "methods" )
Nm <- length( mnam )

dnam <- list( "From:" = mnam,
                "To:" = mnam,
                        c("alpha","beta","sd.pred") )
conv.array <- array( NA, dim=sapply( dnam, length ), dimnames=dnam )
qnt <- t( apply( as.matrix( object ),
                 2,
                 quantile,
                 probs=c(0.5,alpha/2,1-alpha/2) ) )
gtx <- t( apply( as.matrix( object ),
                 2,
                 function(x) c(">0"=mean(x>0,na.rm=TRUE),
                               ">1"=mean(x>1,na.rm=TRUE)) ) )
medians <- qnt[,"50%"]
# Construct the conversion array:
for( ff in 1:Nm ) for( tt in 1:Nm )
   if( ff != tt )
     {
     conv.array[ff,tt,] <-
     medians[c(paste(  "alpha[",mnam[tt],".",mnam[ff],"]",sep=""),
               paste(   "beta[",mnam[tt],".",mnam[ff],"]",sep=""),
               paste("sd.pred[",mnam[tt],".",mnam[ff],"]",sep=""))]

     }
     else
     conv.array[ff,tt,] <-
     c( 0, 1,
     medians[paste("sd.pred[",mnam[tt],".",mnam[ff],"]",sep="")] )
   
# Correction of the median intercepts to make translation
# formulae that are the same both ways (i.e combine to the identity):
for( ff in 1:Nm ) for( tt in 1:Nm )
   if( ff != tt )
     {
     aft <- conv.array[ff,tt,1]
     bft <- conv.array[ff,tt,2]
     atf <- conv.array[tt,ff,1]
     conv.array[ff,tt,1] <-  (  aft      -atf*bft )/2
     conv.array[tt,ff,1] <-  ( -aft/bft + atf     )/2
     }
     
# The variance components
wh <- grep( "sigma", rownames( qnt ) )
wh <- wh[order( rownames( qnt )[wh] )]
var.comp <- qnt[wh,]

# The mean value parameters
alphas <- cbind( qnt[grep("alpha",rownames(qnt)),],
                 gtx[grep("alpha",rownames(qnt)),">0"] )
 betas <- cbind( qnt[grep( "beta",rownames(qnt)),],
                 gtx[grep( "beta",rownames(qnt)),">1"] )
colnames( alphas )[4] <-
colnames(  betas )[4] <- "P(>0/1)"
mean.par <- rbind( alphas, betas )

invisible( list( conv.array = conv.array,
                   var.comp = var.comp,
                   mean.par = mean.par ) )
}

print.MethComp <-
function( x,
     across,
     digits = 3,
      alpha = 0.05,
        ... )
{
Got.coda <- require( coda )
if( !Got.coda )
  stop( "Using the print.MethComp function requires that\n",
        "the package 'coda' is installed.\n",
        "All installed packages are shown if you type 'library()'." )

# Check
if( !inherits( x, "MethComp" ) )
    stop( "\nThe argument to print.MethComp must be of class MethComp." )
# How should we print the conversion table
if( missing( across ) ) across <- ( length( attr( x, "methods" ) ) < 3 )

# Compute the summary
summx <- summary.MethComp( x, alpha=alpha, ... )

# Print a nice summary of the conversion formulae
cat( "\nConversion formula:\n y_to = alpha + beta * y_from +/- 2*sd.pred:\n\n" )
print( round( ftable( summx$conv.array,
                      row.vars = if(across) 2 else c(2,3) ),
                      digits ) )
cat( "\nVariance components with", (1-alpha)*100, "% cred.int.:\n" )
print( round( summx$var.comp, digits ) )
cat( "\nMean parameters with", (1-alpha)*100, "% cred.int.:\n" )
print( round( summx$mean.par, digits ) )
cat("\n Note that intercepts in conversion formulae are adjusted to get",
    "\n conversion fromulae that represent the same line both ways,",
    "\n - therefore are the median of the alphas above not identical",
    "\n to the intercepts given in the conversion formulae.\n")
}
