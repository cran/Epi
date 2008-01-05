apc.plot <-
function( obj, r.txt="Rate", ... )
{
if( !inherits( obj, "apc" ) ) stop( "Argument must be an apc-object" )

# Determine the ranges of the horizontal axes
 a.lab = nice( obj$Age[,1] )
cp.lab = nice( c(obj$Per[,1],obj$Coh[,1]) )[-1]
# The necessary range of the two vertical axes
 r.rg  <- range( obj$Age[,-1] )
rr.rg  <- range( rbind( obj$Per[,-1], obj$Coh[,-1] ) )
# Align the RR with the rates on an integer power of 10
rr.ref <- 10^floor( log10(r.rg[2])-log10(rr.rg[2]) )
# Find the tic-marks for the two vertical axes
 r.tic <- nice(  r.rg, log=T, lpos=1:9 )
rr.tic <- nice( rr.rg, log=T, lpos=1:9 )
# Expand to cover it all
 r.tic <- sort( unique( c( r.tic, rr.tic*rr.ref ) ) )
rr.tic <- r.tic/rr.ref
# Find the places for labels
 r.lab <- nice(  r.tic, log=T, lpos=c(1,2,5) )
rr.lab <- nice( rr.tic, log=T, lpos=c(1,2,5) )
 r.lab <-  r.lab[ r.lab>min( r.tic) &  r.lab<max( r.tic)]
rr.lab <- rr.lab[rr.lab>min(rr.tic) & rr.lab<max(rr.tic)]
# Now draw the frame
fpar <-
apc.frame( a.lab=a.lab,
          cp.lab=cp.lab,
           r.lab=r.lab,
           r.tic=r.tic,
          rr.lab=rr.lab,
          rr.tic=rr.tic,
          rr.ref=rr.ref,
           r.txt=r.txt )
# - and the reference line           
segments( min(cp.lab)-fpar[1], fpar[2], max(cp.lab)-fpar[1], fpar[2] )
apc.lines( obj, frame.par=fpar, ... )
fpar
}
