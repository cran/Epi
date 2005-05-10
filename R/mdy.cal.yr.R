mdy.cal.yr <-
  function( m, d, y )
  {
    if ( is.numeric( m ) )
       {
       if( !is.numeric( d ) || !is.numeric( y ) )
       stop( "All three arguments must be numeric,\n",
             "or the first must be a character vector.\n" )
       res <- cal.yr( ISOdate( y, m, d ) )
       }
    if ( is.character( m ) )
       {
       dl1 <- substr( d[1], 3, 3 )
       dl2 <- substr( d[1], 6, 6 )
       res <- cal.yr( strptime( m,
                      format=paste( "%m", dl1, "%d", dl1, "%Y",  sep="" ) ) )
       }
    res
  }
