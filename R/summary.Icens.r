summary.Icens <-
function( object, scale=1, ... )
{
  if( attr( object, "model" ) == "MRR" )
    {
      if( is.null( object$rates ) )
        {
          class( object ) <- "glm"
          emat <- ci.lin( object )
        }
      else
        {
          rate.est <- ci.lin( object$rates )
          rate.est[,-(3:4)] <- rate.est[,-(3:4)] * scale
          emat <- rbind( cbind( rate.est, RR=NA )[,c(1:4,7,5:6)],
                        ci.lin( object$cov, Exp=T ) )
        }
    }
  if( attr( object, "model" ) == "AER" )
    {
      rate.est <- ci.lin( object$rates )
      rate.est[,-(3:4)] <- rate.est[,-(3:4)] * scale
      emat <- rate.est
    }
  if( length( object ) == 4 )
    {
      b.est <- object[["boot.ci"]]
      colnames( b.est ) <- c( "Boot-med", "lo", "hi" )
      emat <- cbind( emat, b.est )
    }
  emat
}
