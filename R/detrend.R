thinCol <-
function( A, tol = 1e-06, col.num=FALSE )
{
# Remove linearly dependent columns from a matrix
QR <- qr(A, tol = tol, LAPACK = FALSE)
# what columns are linearly dependen on the previous
whCol <- QR$pivot[seq(length = QR$rank)]
# return only numbers of dependent columns if required
if( col.num ) return( whCol )                 
A[, whCol, drop = FALSE]
}

in.span <-
inSpan <-
function( A, x, coef=FALSE, tol=1e-8 )
{
if( is.vector(x) ) dim(x) <- c(length(x),1)
if( nrow(A)!=nrow(x) ) stop("Matrices must have same row dimension")    
# Check if x is in span(A) using regression
xcf  <- NULL
insp <- TRUE
for( i in 1:ncol(x) )
   { 
   mod  <- lm( x[,i] ~ A - 1 )
   insp <- insp & ( all( abs(mod$residuals)<tol ) )
   xcf  <- cbind( xcf, coef(mod) )
   }
xcf[is.na(xcf)] <- 0
if( coef & insp ) {
       return( round( xcf, floor(abs(log10(tol))) ) )
} else return( insp )
}

id.span <-
idSpan <-
function( A, B, tol=1e-8 ) in.span(A,B) & in.span(B,A)
    
detrend <-
function( M, t, weight=rep(1,nrow(M)) )
{
# Detrend the matrix using a weighted inner product.
thinCol( projection.ip( cbind( 1, t ), M , orth = TRUE, weight = weight ) )
}

decurve <-
function( M, t, weight=rep(1,nrow(M)) )
{
# De-trend and -curve the matrix using a weighted inner product.
thinCol( projection.ip( cbind( 1, t, t^2 ), M , orth = TRUE, weight = weight ) )
}
