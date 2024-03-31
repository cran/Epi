print.Icens <-
function( x, scale = 1, digits = 4, ... )
{
emat <- summary.Icens( x, scale=scale )
print( round( emat, digits ) )
}
