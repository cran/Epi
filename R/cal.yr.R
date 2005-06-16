cal.yr <-
function (x, format) 
{
# Check if the input is of some kind of date or time object
  if( any( inherits( x,
                     c("Date","POSIXct","POSIXlt","date","dates","chron")
                    ) ) )
           x <- as.Date( as.POSIXct( x ) )
  else if( is.character( x ) ) {
           if( missing( format ) ) stop( "Character input requires a format." )
           x <- as.Date( x, format = format ) }
  else if( is.factor( x ) ) {
           if( missing( format ) ) stop( "Factor input requires a format." )
           x <- as.Date( as.character( x ), format = format ) }
  else stop( "\nInput should be either character, factor or",
             "some kind of date or time object:\n",
             "Date, POSIXct, POSIXlt, date, dates or chron" )
  res <- as.numeric( x ) / 365.25 + 1970
  res
}
