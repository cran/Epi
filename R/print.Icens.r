print.Icens <-
function( x, digits=4, scale=1, ... )
{
if( attr( x, "model" ) == "MRR" )
{
if( is.null( x$rates ) )
  {
  class( x ) <- "glm"
  rate.est <- ci.lin( x )
  print( round( rate.est, digits ) )
  }
else
  {
  rate.est <- ci.lin( x$rates )
  rate.est[,-(3:4)] <- rate.est[,-(3:4)] * scale
  print( round( rbind( cbind( rate.est, RR=NA )[,c(1:4,7,5:6)],
                       ci.lin( x$cov, E=T ) ),
                digits ) )
  }
}
if( attr( x, "model" ) == "AER" )
{
class( x ) <- "glm"
rate.est <- ci.lin( x )
rate.est[,-(3:4)] <- rate.est[,-(3:4)] * scale
print( round( rate.est, digits ) )
}
}
