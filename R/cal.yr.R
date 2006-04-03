cal.yr <-
function (x, format="%Y-%m-%d") 
{
# Check if the input is some kind of date or time object
  if( any( inherits( x,
                     c("Date","POSIXct","POSIXlt","date","dates","chron")
                    ) ) )
           x <- as.Date( as.POSIXct( x ) )
  else if( is.character( x ) ) x <- as.Date( x, format = format )
  else if( is.factor( x ) ) x <- as.Date( as.character( x ), format = format )
  else stop( "\nInput should be either character, factor or",
             "some kind of date or time object:\n",
             "Date, POSIXct, POSIXlt, date, dates or chron" )
  res <- as.numeric( x ) / 365.25 + 1970
  class( res ) <- c("cal.yr","numeric")
  res
}
