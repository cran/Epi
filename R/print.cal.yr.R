print.cal.yr <-
  function( x, format="%d/%m/%Y", ... )
{
  spyr <- 365 * 24 * 60 * 60
  x <- ( x - 1970 ) * spyr
  class( x ) <- "POSIXct"
  print( strftime( as.POSIXlt( x ), format=format ) )
  invisible( )
}
