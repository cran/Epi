dmy.cal.yr <-
  function( d, m, y )
  {
    if ( is.numeric( d ) )
       {
       if( !is.numeric( m ) || !is.numeric( y ) ) {
       stop( "All three arguments must be numeric,\n",
             "or the first must be a character vector.\n" ) }
       res <- cal.yr( ISOdate( y, m, d ) )
       }
 if ( is.character( d ) )
       {
       dl1 <- substr( d[1], 3, 3 )
       dl2 <- substr( d[1], 6, 6 )
       res <- cal.yr( strptime( d,
                      format=paste( "%d", dl1, "%m", dl1, "%Y",  sep="" ) ) )
       }
    res
  }
