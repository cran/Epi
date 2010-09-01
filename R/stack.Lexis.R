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
function( x, ... )
{
zz <- table(x$lex.Cst,x$lex.Xst)
class(zz) <- "matrix"
diag(zz) <- NA
zz[zz==0] <- NA
zz
}

# The factorize method
factorize <- function (obj, ...) UseMethod("factorize")

factorize.Lexis <-
function( obj, ... )
{
   obj$lex.Cst <- factor( obj$lex.Cst )
   obj$lex.Xst <- factor( obj$lex.Xst )
   all.levels = union(levels(obj$lex.Cst),levels(obj$lex.Xst))
   obj$lex.Cst <- factor( obj$lex.Cst, levels=all.levels )
   obj$lex.Xst <- factor( obj$lex.Xst, levels=all.levels )
   obj
}

# The msdata method
msdata <- function (obj, ...) UseMethod("msdata")

msdata.Lexis <-
function( obj,
   time.scale = timeScales(obj)[1], ... )
{
if( !require( mstate ) )
  stop( "You do not want this before you have installed the 'mstate' package.\n" )
tr.mat <- tmat(obj)
# Essentially a msdata object is a stacked Lexis object with other variable names
tmp <- stack.Lexis( factorize.Lexis(obj) )
lv  <- c( match(timeScales(obj), names(tmp) ),
          grep("lex\\.", names(tmp) ) )
# The resulting dataframe is created by renaming columns in the stacked Lexis object
data.frame( id = tmp$lex.id,
          from = as.integer( tmp$lex.Cst ),
            to = as.integer( tmp$lex.Xst ),
         trans = as.integer( tmp$lex.Tr ),
        Tstart = tmp[,time.scale],
         Tstop = tmp[,time.scale] + tmp$lex.dur,
          time = tmp$lex.dur,
        status = as.integer( tmp$lex.Fail ),
                 tmp[,-lv] )
}
