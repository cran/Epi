cal.yr <-
function( x,
     format = "%Y-%m-%d",
         wh = NULL )
{
cl.typ <- c("Date","POSIXct","POSIXlt","date","dates","chron")
# Check if the input is a data frame and convert column bu column
if (inherits(x, "data.frame"))
   { 
   if (is.null(wh) & missing(format))
      {
      # Indicator of where a date-type variable is
      wh <- which(sapply(x, inherits, cl.typ))
      }
   if (is.null(wh) & !missing(format))
      {
      # Indicator of where the character variables are
      wh <- which(sapply(x, is.character))
      }
   if (is.character(wh)) wh <- match(wh, names(x))
   # Convert the dates or the character variables one at a time
   for (i in wh)
       {
       if (is.character(x[,i]))
          x[,i] <- cal.yr(x[,i], format = format)
       else
          x[,i] <- cal.yr(x[,i])
       }
   return(x)
   } # end of dataframe doings
# now converting a single vector to decimal years:
# Check if the input is some kind of date or time object
if (any(inherits(x, cl.typ))) x <- as.Date(as.POSIXct(x))
else if( is.character(x))     x <- as.Date(x,
                                           format = format)
else if( is.factor(x))        x <- as.Date(as.character(x),
                                           format = format)
else stop("\nInput should be a data frame, a character vector",
          "a factor or some kind of date or time object:\n",
          "Date, POSIXct, POSIXlt, date, dates or chron")
res <- as.numeric(x) / 365.25 + 1970
class(res) <- c("cal.yr", "numeric")
return(res)
}
