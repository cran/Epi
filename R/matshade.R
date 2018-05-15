matshade <-
function( x, y,
          lty = 1,
          col = 1:(ncol(y)/3), col.shade=col, alpha=0.15,
         plot = dev.cur()==1,
          ... )
    {
# check sanity of x and y
if( !is.vector(x) ) stop( "x must be a vector\n")        
if( ncol(y)%%3 != 0 ) warning( "number of columns in y, ", ncol(y),
                               " is not divisible by 3\n")
# create a new plot?
if( plot ) matplot( x, y, col="transparent", ... )
# number of curves to draw        
ncrv <- ncol(y) %/% 3
# The recycling rule for colors and alpha:
col       <- rep(col      ,ncrv)[1:ncrv]
col.shade <- rep(col.shade,ncrv)[1:ncrv]
alpha     <- rep(alpha    ,ncrv)[1:ncrv]
# First shaded areas        
for( i in 1:ncrv ) polygon( c( x        ,rev(x)        ),
                            c( y[,i*3-1],rev(y[,i*3-0])),
                           col = adjustcolor( col.shade[i],
                                          alpha.f=alpha[i] ),
                           pch = NULL,
                        border = "transparent" )
# then curves on top of these
matlines( x, y[,(1:ncrv)*3-2], col=col, lty=lty, ... )
    }
