# The method
nid <- function ( Lx, ... ) UseMethod("nid")

nid.default <-
nid.Lexis <-
function( Lx, by=NULL, ... )
{
if( !is.null(by) )
  {
  if( !(by %in% names(Lx)) ) stop( "'by' must be the name of a variable in Lx" )
  if( !is.factor(Lx[,by]) ) Lx[,by] <- factor(Lx[,by])
  }
if( is.null(by) )
  length( unique( Lx$lex.id ) )
else
  sapply( split( Lx, Lx[,by] ), nid.Lexis ) 
}

# Make a boostrap sample of a Lexis object:
# Sample the *persons* with replacement, possibly sampling within levels of by=
bootLexis <-
function( Lx,
        size = NULL,
          by = NULL,
     replace = TRUE )
{
if( !inherits( Lx, "Lexis" ) ) stop("Only meaningful for Lexis objects.")
    
isDT <- inherits( Lx, "data.table" )
if( isDT ) class( Lx ) <- c("Lexis","data.frame")
    
# determine size of the bootstrap samples if not given   
if( is.null( size ) ) size <- nid.Lexis( Lx, by = by )

# allowing for a length 1 x-vector
REsample <- function(x,sz) x[sample.int(length(x),size=sz,replace=replace)]
    
if( is.null(by) ) { # Simple bootstrap
  bLx <- subid.Lexis( Lx, REsample( unique(Lx$lex.id), size ) )
} else { # Bootstrap by groups
  bLx <- NULL
  spL <- split( Lx, Lx[,by] )
  for( i in 1:length(spL) )
     {
     bLx <- rbind( bLx,
                   cbind( bootLexis( spL[[i]], size = size[i] ),
                          bgr = paste(i) ) )
     }
  bLx$lex.id <- as.integer( interaction(bLx$lex.id,bLx$bgr) )
  bLx <- bLx[,-grep("bgr",names(bLx))]
  }
# return the result after converting to data.table if needed
if( isDT ) class( bLx ) <- c("Lexis","data.table","data.frame")
bLx
}

# A utility function that returns a Lexis object subsetted to a set of
# lex.ids, allowing for repeat values of lex.id
subid.Lexis <-
function( Lx, ids )
{
tt <- table( ids )
bLx <- NULL
max.id <- 0
for( j in 1:max(tt) )
   {
   # avoid note about no visible binding
       wh <- NULL
   lex.id <- NULL
   # who appears at least j times in the sample ?
   wh <<- names(tt[tt>=j])
   sb <- subset( Lx, lex.id %in% wh )
   # remember their original id
   sb$old.id <- sb$lex.id
   # to assure that different samples of the same person has different lex.id
   sb$lex.id <- sb$lex.id + max.id
   max.id <- max( sb$lex.id )
   bLx <- rbind( bLx, sb )
   }   
# Generate new lex.ids in the order 1:N
bLx$lex.id <- as.integer( factor(bLx$lex.id) )
bLx
}
