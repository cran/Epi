# Functions to facilitate analysis of multistate models

# The stack method is already defined (in the utils package)
# and hence imported in the NAMESPACE file

stack.Lexis <-
function( x, ... )
{
## Function to stack obervations for survival analysis
## Same covariates
xx <- data.frame( cbind( x, lex.Tr="", lex.Fail=FALSE ) )[NULL,]
tm <- tmat.Lexis( x )
for( fst in levels(factor(x$lex.Cst)) )
for( tst in levels(factor(x$lex.Xst)) )
   if( !is.na(tm[fst,tst]) )
     {
     tr = paste( fst, "->", tst , sep="" )
     zz <- x[x$lex.Cst==fst,]
     xx <- rbind( xx, data.frame( zz, lex.Tr=tr, lex.Fail=(zz$lex.Xst==tst) ) )
     }
xx$lex.Tr <- factor(xx$lex.Tr)
## Reshuffle the variables
wh.om <- match( "lex.Xst", names(xx) )
wh.rl <- match( c("lex.Tr","lex.Fail"), names(xx) )
xx <- xx[,c(1:wh.om,wh.rl,(wh.om+1):(wh.rl[1]-1))]
class( xx ) <- c("stacked.Lexis","data.frame")
xx
}

# The tmat method
tmat <- function (x, ...) UseMethod("tmat")

tmat.Lexis <-
function( x, Y=FALSE, mode="numeric", ... )
{
zz <- table(x$lex.Cst,x$lex.Xst)
class(zz) <- "matrix"
if( Y )
  {
  diag(zz) <- summary( x, simplify=FALSE )[[1]][1:nrow(zz),"Risk time:"]
  }
else diag(zz) <- NA
zz[zz==0] <- NA
if( mode != "numeric" ) zz <- !is.na(zz)
zz
}
