cal.yr <-
function (x, format) 
{
    if (mode(x) == "character") x <- strptime( x, format = format)
    res <- as.numeric( as.Date( x ) ) / 365.25 + 1970
    attributes( res ) <- NULL
    attributes( res )$class <- "cal.yr"
    res
}
