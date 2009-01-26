BA.plot <-
function( y1, y2,
      meth.names = NULL,
       mean.repl = FALSE,
     comp.levels = 1:2,
             ... )
{
  if( is.data.frame( y1 ) )
    {
    # If in the long form, convert to 2-column matrix
    if( sum( !is.na( match( c("item","meth","y"), names( y1 ) ) ) ) == 3  )
      {
      repl <- !is.na( match( c("repl"), names( y1 ) ) )
      # Remove any disturbing empty factor levels
      if( is.null( meth.names ) ) meth.names <- levels( factor( y1$meth ) )
      y1$meth <- as.factor( as.integer( y1$meth ) )
      yy <- as.data.frame( tapply( y1[,"y"],
                                   c( if( repl & !mean.repl )
                                      list( interaction( y1[,"item"],
                                                         y1[,"repl"] ) )
                                      else list( y1[,"item"] ),
                                      list( y1[,"meth"] ) ),
                                   mean ) )
      yy <- yy[,comp.levels]
      yy <- yy[complete.cases(yy),]
      if( nrow(yy)==0 ) stop( "No items have measurements by both methods '",
                              meth.names[comp.levels[1]], "' and '",
                              meth.names[comp.levels[2]], "'." )
      BA.plot( yy, meth.names=meth.names[comp.levels], ... )
      }
    else
    # If a two-column matrix
      {
      if( dim(y1)[2]==2 )
        {
        meth.names <- if( is.null( meth.names ) ) names( y1 )
                      else meth.names
        y2 <- y1[,2]
        y1 <- y1[,1]
        BA.plot( y1, y2, meth.names=meth.names, ... )
        }
      }
    }
  else
  # If two vectors are supplied
    {
    if( is.null( meth.names ) )
      {
      n1 <- deparse( substitute( y1 ) )
      n2 <- deparse( substitute( y2 ) )
      }
    else
      {
      n1 <- meth.names[1]
      n2 <- meth.names[2]
      }
    BlandAltman( y1, y2, x.name=n1, y.name=n2, ... )
    }
}
