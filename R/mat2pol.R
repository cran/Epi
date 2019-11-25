mat2pol <-
function( pm,
        perm = 1:ncol(pm),
           x = as.numeric(rownames(pm)),
         col = rainbow(ncol(pm)),
          yl = 0:1,
      append = FALSE,
         ... )
{
if( missing(x) ) x <- as.numeric( rownames(pm) )
if( length(x)==0 ) x <- 1:nrow(pm)
xm <- cbind( 0, pm[,perm] )
xm <- t( apply(xm,1,cumsum) )
if( !append ) plot( x, x, type="n", ylim=yl, ... )
for( j in 1:ncol(pm) )
   polygon( c(x,rev(x)),
            c(xm[,j],rev(xm[,j+1])),
            col = col[j], border="transparent" )
invisible( xm )
}
