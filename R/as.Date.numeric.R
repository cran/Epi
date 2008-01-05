as.Date.numeric <-
function( x, ..., unit="d" )
{
if(!(toupper( substr( unit, 1, 1 ) ) %in% c("Y","D") ) )
  stop( "2nd argument must be either 'd' or 'y'" )
if( toupper( substr( unit, 1, 1 ) ) == "Y" | inherits( x, "cal.yr" ) )
  xx <- structure( round( ( x - 1970 ) * 365.25 ), class="Date" )
if( toupper( substr( unit, 1, 1 ) ) == "D" )
  xx <- structure( x, class="Date" )
xx
}
