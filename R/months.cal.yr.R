months.cal.yr <-
function( x, abbreviate=FALSE )
{
spyr <- # Seconds per year
  difftime( strptime( "01/01/1972", format="%d/%m/%Y" ),
            strptime( "01/01/1968", format="%d/%m/%Y" ), units="secs"
           ) / 4
mon <- rep( strptime( "01/01/1952", format="%d/%m/%Y" ), 12 )
mon$mon <- mon$mon + 0:11
mnam <- months( mon, abbreviate=abbreviate )
z <- ( x - 1970 ) * spyr
class( z ) <- c( "POSIXt", "POSIXct" )
nam <- format(z, ifelse(abbreviate, "%b", "%B") )
factor( nam, levels=mnam )
}
