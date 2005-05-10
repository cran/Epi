weekdays.cal.yr <-
  function( x, abbreviate=FALSE )
{
spyr <- # Seconds per year
  difftime( strptime( "01/01/1972", format="%d/%m/%Y" ),
            strptime( "01/01/1968", format="%d/%m/%Y" ), units="secs"
           ) / 4
# Get the weekday names to use as levels of the factor
sun <- rep( strptime( "13/07/1952", format="%d/%m/%Y" ), 7 )
sun$wday <- sun$wday + 0:6
wnam <- weekdays( sun, abbreviate=abbreviate )
# Then compute the weekdays and attach the proper factor levels
z <- ( x - 1970 ) * spyr
class( z ) <- c( "POSIXt", "POSIXct" )
nam <- format(z, ifelse( abbreviate, "%a", "%A") )
factor( nam, levels=wnam )
}
