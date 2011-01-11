# The factorize method
factorize <- function (obj, ...) UseMethod("factorize")

# The Lexis version of this
factorize.Lexis <-
function (obj, states=NULL, print=TRUE, ...)
{
# Is this really a Lexis object
if( !inherits(obj,"Lexis") ) stop( "First argument must be a Lexis object" )

# Make the state variables into factors with the same levels set
Cst <- factor(obj$lex.Cst)
Xst <- factor(obj$lex.Xst)
all.levels = union(levels(Cst), levels(Xst))
Cst <- factor(Cst, levels = all.levels)
Xst <- factor(Xst, levels = all.levels)
obj$lex.Cst <- Cst
obj$lex.Xst <- Xst

# If new state names as a list it implies merging of them
if( !is.null( states ) & is.list( states ) )
  {
  # Collapsing levesl requires the input factor to have levels that
  # are given in the elements of the list
  levels( obj$lex.Cst ) <- unlist(states)
  levels( obj$lex.Xst ) <- unlist(states)
  obj$lex.Cst <- Relevel( obj$lex.Cst, states, ... )
  obj$lex.Xst <- Relevel( obj$lex.Xst, states, ... )
  if( print )
    {
    # Construct translation table between old and grouped states to print
    tC <- table( Cst, obj$lex.Cst )
    tX <- table( Xst, obj$lex.Xst )
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
  return( obj )
  }

# If states are just given as a vector we assume that it's just new names
if( !is.null( states ) & !is.list( states ) )
  {
  if( length( states ) != nlevels(Cst) )
    stop( "Second argument is a vector of length ", length(states),
          ", but it should be the joint no. of states, ",
         length(all.levels),
          "\ncorresponding to ", all.levels )
  levels( obj$lex.Cst ) <-
  levels( obj$lex.Xst ) <- states
  if( print )
    {
    cat( "New levels for lex.Xst and lex.Cst generated:\n" )
    print( data.frame( old=all.levels, new=levels(obj$lex.Cst) ) )
    }
  }

return( obj )
}
