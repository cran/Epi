# The factorize method
factorize <- function (x, ...) UseMethod("factorize")

# Default method is just the Relevel method
factorize.default <- Relevel.default

# The Lexis version of this
  Relevel.Lexis <-
factorize.Lexis <-
function (x, states=NULL, print=TRUE, ...)
{
# Is this really a Lexis object
if( !inherits(x,"Lexis") ) stop( "First argument must be a Lexis object" )

# Make the state variables into factors with the same levels set
Cst <- factor(x$lex.Cst)
Xst <- factor(x$lex.Xst)
all.levels = union(levels(Cst), levels(Xst))
Cst <- factor(Cst, levels = all.levels)
Xst <- factor(Xst, levels = all.levels)
x$lex.Cst <- Cst
x$lex.Xst <- Xst

# If new state names are given as a list it implies merging of them
if( !is.null( states ) & is.list( states ) )
  {
  x$lex.Cst <- Relevel( x$lex.Cst, states, ... )
  x$lex.Xst <- Relevel( x$lex.Xst, states, ... )
  if( print )
    {
    # Construct translation table between old and grouped states to print
    tC <- table( Cst, x$lex.Cst )
    tX <- table( Xst, x$lex.Xst )
    cC <- matrix( colnames(tC), nrow(tC), ncol(tC), byrow=T )
    cX <- matrix( colnames(tX), nrow(tX), ncol(tX), byrow=T )
    cC[tC==0] <- ""
    cX[tX==0] <- ""
    print( data.frame( type=rep( c("lex.Cst","lex.Xst"),
                                 c(nrow(tC),nrow(tX)) ),
                       old=c(rownames(tC),rownames(tX)),
                       new=c( apply( cC, 1, paste, collapse="" ),
                              apply( cX, 1, paste, collapse="" ) ) ) )
    }
  return( x )
  }

# If states are just given as a vector we assume that it's just new names
if( !is.null( states ) & !is.list( states ) )
  {
  if( length( states ) != nlevels(Cst) )
    stop( "Second argument is a vector of length ", length(states),
          ", but it should be the joint no. of states, ",
         length(all.levels),
          "\ncorresponding to ", all.levels )
  levels( x$lex.Cst ) <-
  levels( x$lex.Xst ) <- states
  if( print )
    {
    cat( "New levels for lex.Xst and lex.Cst generated:\n" )
    print( data.frame( old=all.levels, new=levels(x$lex.Cst) ) )
    }
  }

return( x )
}
