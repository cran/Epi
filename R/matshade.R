matshade <-
function( x, y,
          lty = 1,
          col = 1:(ncol(y)/3),
    col.shade = col,
        alpha = 0.15,
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
# First shaded areas - allowing NA to create holes
draw.shade <-
function( xx, yu, yl, col )
{
# Find NAs and plot shades for non-missing sequences
wh <- 1:length(xx)
# NAs where at least one is missing
wh[apply( cbind(xx,yl,yu), 1, function(x) any(is.na(x)) )] <- NA
idx <- 1 + cumsum( is.na( wh ) )
not.na <- !is.na( wh )
chunks <- split( wh[not.na], idx[not.na] )
for( ch in chunks )
if( length(ch) > 1 ) polygon( c( xx[ch], rev(xx[ch]) ),    
                              c( yu[ch], rev(yl[ch]) ),    
                              col = col,           
                              pch = NULL,          
                           border = "transparent" )
}
# Then the shades for each set of 3 columns
for( i in 1:ncrv ) draw.shade( x, y[,i*3-1],
                                  y[,i*3-0], 
                               col = adjustcolor( col.shade[i],
                                              alpha.f=alpha[i] ) )
# then curves on top of these
matlines( x, y[,(1:ncrv)*3-2], col=col, lty=lty, ... )
    }
       
