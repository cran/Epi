weeks.cal.yr <-
  function( x )
{
spyr <- # Seconds per year
  difftime( strptime( "01/01/1972", format="%d/%m/%Y" ),
            strptime( "01/01/1968", format="%d/%m/%Y" ), units="secs"
           ) / 4
z <- ( x - 1970 ) * spyr
class( z ) <- c( "POSIXt", "POSIXct" )
wkno <- format(z, "%W" )
as.numeric( wkno )
}
