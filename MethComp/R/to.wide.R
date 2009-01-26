to.wide <-
function( data )
{
if( !is.data.frame(data) )
stop( "The argument must be a dataframe\n--- you supplied a ", class(data) )

# Check names of the supplied dataframe
rq.nam <- c("meth", "item", "y")
if (sum(!is.na(wh <- match(rq.nam, names(data)))) < 3)
    stop("\nThe supplied dataframe misses column(s) named ",
          rq.nam[is.na(wh)], ".\n")

# Are replicates present ?
if( "repl" %in% names( data ) )
  {
  data$id <- interaction( data$item, data$repl )
  if( length( unique( data$id ) ) > length( unique( data$item ) ) )
  cat( "\nNote:\n Replicate measurements are taken as separate items!\n" )
  }
else
  data$id <- data$item

# Method must be a factor with only the levels actually present
data$meth <- factor( data$meth )

res <- reshape( data, direction = "wide",
                        v.names = "y",
                        timevar = "meth",
                          idvar = "id" )
names( res ) <- gsub( "y\\.", "", names( res ) )
attr( res, "reshapeWide" )$varying <- gsub( "y\\.", "",
                                            attr( res, "reshapeWide" )$varying )
res
}

to.long <-
function( data, vars )
{
if( !is.data.frame(data) )
  stop( "The argument must be a dataframe\n--- you supplied a ", class(data) )
if(      missing( vars ) )
  stop( "You must supply names or numbers of the variables to be stacked" )
if( is.numeric  ( vars ) ) v.names <- names( data )[vars]
if( is.character( vars ) ) v.names <- vars
mn <- match( v.names, names( data ) )
if( any( is.na( mn ) ) ) stop( v.names[is.na( mn )], " is not in the dataframe" )
if( missing( vars ) ) stop("You must specify which variables contain measurements")
new <-
reshape( data, direction = "long",
                 varying = list(v.names),
                   times = v.names,
                 v.names = "y",
                 timevar = "meth",
                   idvar = "item" )
rownames( new )  <- NULL
new$meth <- factor( new$meth, levels=v.names )
new
}
