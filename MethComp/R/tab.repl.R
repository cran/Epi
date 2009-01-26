tab.repl <-
function( data )
{
if( !is.data.frame( data ) )
  stop("\nThe argument must be a dataframe, not a ", class( data ) )
rq.nam <- c("meth","item","repl","y")
if( sum( !is.na( wh <- match( rq.nam, names( data ) ) ) ) < 4 )
  stop("\nThe supplied dataframe misses column(s) named ",
       rq.nam[is.na(wh)], ".\n" )
# Make a table of no. replicates for each item
# Table of item by method
tt <- with( data, tapply( !is.na(y), list(item,meth), sum ) )
# What sort of numbers of replicates does actually exist?
xx <- as.numeric( names( table( tt ) ) )
# Add a row of each to tt so that the result of tt is rectangular
tt <- rbind( tt, matrix(xx,length(xx),ncol(tt)) )
# Subtract 1 from the result t compensate for the added rows
XX <- t(apply(tt,2,table)-1)
if( dim(XX)[1]==1 )
  {
  XX <- t(XX)
  colnames(XX) <- xx
  }
# Total no. observations by method
tt <- with( data, tapply( !is.na(y), list(meth), sum ) )
if( dim(XX)[2]>1 ) XX <- addmargins(XX,2)
XX <- cbind( XX, tt )
# Niceify the result
colnames( XX )[ncol(XX)-(1:0)] <- c("#Items",paste("#Measurements:",sum(!is.na(data$y))))
names(dimnames(XX)) <- c("Method"," #Replicates")
XX
}
