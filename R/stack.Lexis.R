# Functions to facilitate analysis of multistate models

# The stack method
stack <- function (x, ...) UseMethod("stack")

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
function( x )
{
zz <- table(x$lex.Cst,x$lex.Xst)
class(zz) <- "matrix"
diag(zz) <- NA
zz[zz==0] <- NA
zz
}

# The mstate method
mstate <- function (obj, ...) UseMethod("mstate")

mstate.Lexis <-
function( obj,
   time.scale = timeScales(obj)[1] )
{
if( !require( mstate ) )
  stop( "You do not want this befor you have installed the 'mstate' package.\n" )
tr.mat <- tmat(obj)
new <- NULL
for( i in 1:nrow(tr.mat) )
for( j in 1:ncol(tr.mat) )
{
if( !is.na(tr.mat[i,j]) )
  {
  # Any transition out of state i generates a record
  tmp <-  obj[obj$lex.Cst==rownames(tr.mat)[i],]
  status <- ( tmp$lex.Xst==colnames(tr.mat)[j] )
  from <- rownames(tr.mat)[i]
  to   <- colnames(tr.mat)[j]
  tmp  <- data.frame( from, to, status, tmp )
  new  <- rbind( new, tmp )
  }
}
id     <-          new[,"lex.id"]
Tstart <-          new[,time.scale]
Tstop  <- Tstart + new[,"lex.dur"]
transition <- interaction( new$from, new$to )
# Dump the old Lexis variables
rm <- c( grep("lex.", names(new) ),
         match( timeScales(obj), names(new) ) )
new <- new[,-rm]
new <- data.frame( id, Tstart, Tstop,
                    from = as.integer(new[,1]),
                      to = as.integer(new[,2]),
                   trans = as.integer( interaction( new$from, new$to ) ),
                   new[,-(1:2)] )
new <- new[order(new$id,new$Tstart),]
# expand.covs( new, tmat(obj), covs, append=TRUE )
# class( new ) <- c("mstate","data.frame")
return( new )
}
