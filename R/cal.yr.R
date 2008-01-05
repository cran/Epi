cal.yr <-
function (x, format="%Y-%m-%d") 
{
cl.typ <- c("Date","POSIXct","POSIXlt","date","dates","chron")
# Check if the input is a data frame and convert
  if( inherits( x, "data.frame" ) )
    {
    # Indicator of where a date-type variable is
    wh <- sapply( x, inherits, cl.typ )
    # The positions
    wh <- (1:length(wh))[wh]
    # Convert them
    for( i in wh ) x[,i] <- cal.yr( x[,i] )
    return( x )
    }
# Check if the input is some kind of date or time object
  if( any( inherits( x, cl.typ ) ) )
           x <- as.Date( as.POSIXct( x ) )
  else if( is.character( x ) ) x <- as.Date( x, format = format )
  else if( is.factor( x ) ) x <- as.Date( as.character( x ), format = format )
  else stop( "\nInput should be a data frame, a character vector, a factor or ",
             "some kind of date or time object:\n",
             "Date, POSIXct, POSIXlt, date, dates or chron" )
  res <- as.numeric( x ) / 365.25 + 1970
  class( res ) <- c("cal.yr","numeric")
  return( res )
}
