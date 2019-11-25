harm <-
function( x, ord=1, per=1, verbose=FALSE )
{
vn <- deparse(substitute(x))
if( inherits(x,"Date") )
  {
  x <- julian( x )
  if( missing(verbose) ) verbose <- TRUE
  if( missing(per) ) per <- 365.25
  if( verbose ) cat( "Note: Date variable", vn, "was converted to days.\n" )
  }
if( verbose ) cat("Note: Period taken to be", per, "\n" )
hm <- NULL
for( i in 1:ord )
   hm <- cbind( hm, sin(2*i*pi*x/per), cos(2*i*pi*x/per) )
colnames(hm) <- paste0( rep(c("sin","cos"),ord), rep(1:ord,each=2) )
return( hm )
}
