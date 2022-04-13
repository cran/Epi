# The factorize method
factorize <- function (x, ...) UseMethod("factorize")

# Default method is just the Lexis method
factorize.default <- 
factorize.Lexis <-
function(x, ..., verbose = FALSE)
{
Cst <- Xst <- NULL
# If lex.Cst and lex.Xst are not factors, make sure they are, and
# if they are restrict levels to to actually ocurring
Cst <- factor(x$lex.Cst)
Xst <- factor(x$lex.Xst)
all.levels <- union(levels(Cst),
                    levels(Xst))
# then amend them to have the same set of levels
x$lex.Cst <- factor(Cst, levels = all.levels)
x$lex.Xst <- factor(Xst, levels = all.levels)

# Note if levels changed
if((length(setdiff(all.levels, levels(Cst))) > 0 |
    length(setdiff(all.levels, levels(Xst))) > 0)
   & verbose)
   cat( "NOTE: lex.Cst and lex.Xst now have levels:\n", all.levels, "\n")  
return(x)
}

# The Lexis version of Relevel
Relevel.Lexis <-
function(x, ref, ...)
{
# Is this really a Lexis object
if(!inherits(x, "Lexis")) stop("First argument must be a Lexis object")
    
# Make sure the states are factors with the same levels
if(!is.factor(x$lex.Cst) | !is.factor(x$lex.Xst))
   stop("lex.Cst and lex.Xst must be factors")
if(any(levels(x$lex.Xst) != levels(x$lex.Cst)))
   stop("lex.Cst and lex.Xst must have same levels\n",
        "you may want to factorize.Lexis()\n")    

# Then use the same Relevel on the factors
x$lex.Cst <- Relevel.factor(x$lex.Cst, ref, ...)
x$lex.Xst <- Relevel.factor(x$lex.Xst, ref, ...)
return(x)
}

# The levels method is already defined (in the base package)
# and hence imported in the NAMESPACE file
levels.Lexis <-
function( x )
{
union(base::levels(x$lex.Cst),
      base::levels(x$lex.Xst))
}

## # The Lexis version of levels
## "levels.Lexis<-" <-
## function(x, value)
## {
## base::levels(x$lex.Cst) <- value
## base::levels(x$lex.Xst) <- value
## x
## }
