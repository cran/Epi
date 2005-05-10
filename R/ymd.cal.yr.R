ymd.cal.yr <-
  function( y, m = NA, d = NA )
  {
    if ( is.numeric( y ) )
       {
       if( !is.numeric( m ) || !is.numeric( d ) )
         stop( "All three arguments must be numeric,\n",
               "or the first must be a character vector.\n" )
       res <- cal.yr( ISOdate( y, m, d ) )
       }
    if ( is.character( y ) )
       {
       dl1 <- substr( y[1], 5, 5 )
       dl2 <- substr( y[1], 8, 8 )
       res <- cal.yr( strptime( y,
                      format=paste( "%Y", dl1, "%m", dl2, "%d",  sep="" ) ) )
       }
    res
  }
