month <-
function( x, abbreviate=FALSE, fact=TRUE )
{
mon <- rep( strptime( "01/01/1952", format="%d/%m/%Y" ), 12 )
mon$mon <- mon$mon + 0:11
mths <- months( mon, abbreviate=abbreviate ) 
if( fact ) factor( months( x ), levels=mnam ) else months( x )
}

weekday <-
function( x, abbreviate=FALSE, fact=TRUE )
{
dnam <- weekdays( as.Date( "07/01/1952",
                  format="%d/%m/%Y" ) + 0:6,
                  abbreviate=abbreviate )
days <- weekdays( x, abbreviate=abbreviate ) 
if( fact ) factor( days, levels=dnam ) else days
}
