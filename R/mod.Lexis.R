modLexis <- 
function( Lx,
      nameLx,
     formula,
        from = preceding(Lx,to),
          to = absorbing(Lx),
      paired = FALSE,
        link = "log", scale = 1, verbose = TRUE,
       model,
         ... )
{
# a common wrapper for glm and gam modeling of Lexis FU  
# is this a Lexis object ?
if( !inherits(Lx,"Lexis") ) stop("The first argument must be a Lexis object.\n")
    
# check that events are actual levels of lex.Xst 
if( is.numeric(to) ) to <- levels( Lx$lex.Xst )[to]
wh <- match( to, levels(Lx$lex.Xst) )
if( any(is.na(wh)) ) stop("'to' must be a subset of: '", 
                     paste(levels(Lx$lex.Xst), collapse="','", sep=""), "'\n" )  

# check that from are actual levels of lex.Cst 
if( is.numeric(from) ) from <- levels( Lx$lex.Cst )[from]
wh <- match( from, levels(Lx$lex.Cst) )
if( any(is.na(wh)) ) stop("'from' must be a subset of: '", 
                     paste(levels(Lx$lex.Cst), collapse="','", sep=""), "'\n" )  
Lx <- Lx[Lx$lex.Cst %in% from,]                     

# work out which transitions are modeled
# first a small utility (transition as text)
trt <- function( f, t ) paste( f, "->", t, sep="" )
if( paired )
  {
if( length(from) != length(to) )
    stop("If 'paired' is TRUE, from and to must have same length!\n")
if( any(from==to) )
    stop("If 'paired' is TRUE, entries in from and to must be different within pairs\n")
trnam <- trt( from, to )
  } else {  
tm <- tmat( Lx )[from,to,drop=FALSE]
trnam <- outer( rownames(tm), colnames(tm), trt )[tm>0]
trnam <- trnam[!is.na(trnam)]
  }
# just for formatting the explanatory text
onetr <- length( trnam )==1
trprn <- paste( trnam, collapse=", " )
    
# warn if a potentially silly model is defined
if( any( (ts<-table(sapply( strsplit(trnam,"->"), function(x) x[1] )))>1 ) ) warning(
 "NOTE:\nMultiple transitions *from* state '",names(ts[ts>1]),"' - are you sure?",
 "\nThe analysis requested is effectively merging outcome states.", 
 "\nYou may want analyses using a *stacked* dataset - see ?stack.Lexis\n" )

# Beginning of a new feature with a countmultiplier of the transitions
# allowing tabular records to be merged to a Lexis object
# --- not used subsequently in this function (yet)
# if( !("lex.N" %in% names(Lx)) ) Lx$lex.N <- 1
 
# construct the model formula - note that we already made sure that
# from and to are pairwise different
if( length(formula) != 2 ) stop("formula must be a one-sided formula")

form <- cbind( trt(Lx$lex.Cst,Lx$lex.Xst) %in% trnam, #*Lx$lex.N,
               Lx$lex.dur ) ~ 1
## form <- cbind( (Lx$lex.Xst %in% to &
##                 Lx$lex.Xst != Lx$lex.Cst), #*Lx$lex.N,
##                 Lx$lex.dur ) ~ 1
form[3] <- formula[2]
from <- levels( factor(Lx$lex.Cst) ) # only levels present in lex.Cst

# Scaling
Lx$lex.dur <- Lx$lex.dur/scale
    
# Tell what we intend to and then do it
if( verbose ){
cat( deparse(substitute(model)),
     " Poisson analysis of Lexis object ", nameLx, " with ", link, " link",
     ":\nRates for", if(  onetr ) " the", " transition",
                     if( !onetr ) "s", ": ", trprn,
     if( scale!=1 ) paste(", lex.dur (person-time) scaled by", scale ), "\n", sep="" )
             }
    
# Fit the model
mod <- model( form, family = poisreg(link=link), data = Lx, ... )
    
# Add an explanatory attribute
attr( mod, "Lexis" ) <- list( data=nameLx,
                             trans=trnam,
                           formula=form[-2],
                             scale=scale )     
mod
}

# Here are the actual functions of interest:
######################################################################
# the glm function
glm.Lexis <- 
function( Lx,
     formula,
        from = preceding(Lx,to),
          to = absorbing(Lx),
      paired = FALSE,
        link = "log",
       scale = 1,
     verbose = TRUE,
         ... )
{
# name of the supplied object
nameLx <- deparse(substitute(Lx))

# sensible defaults if one of to and from is missing
if(  missing(from) & !missing(to) ) from <- preceding (Lx,to  )
if( !missing(from) &  missing(to) ) to   <- succeeding(Lx,from)
xx <- modLexis( Lx, nameLx,
                formula, from, to,
                paired = paired, link = link, scale = scale, verbose = verbose,
                 model = stats::glm, ... )
class( xx ) <- c( "glm.lex", class(xx) )
xx
}

######################################################################
# the gam function
gam.Lexis <- 
function( Lx,
     formula,
        from = preceding(Lx,to),
          to = absorbing(Lx),
      paired = FALSE,
        link = "log",
       scale = 1,
     verbose = TRUE,
         ... )
{
# name of the supplied object
nameLx <- deparse(substitute(Lx))

# sensible defaults if one of the two is missing
if(  missing(from) & !missing(to) ) from <- preceding (Lx,to  )
if( !missing(from) &  missing(to) ) to   <- succeeding(Lx,from)
xx <- modLexis( Lx, nameLx,
                formula, from, to,
                paired = paired, link = link, scale = scale, verbose = verbose,
                 model = mgcv::gam, ... )
class( xx ) <- c( "gam.lex", class(xx) )
xx
}

######################################################################
# And here is the coxph counterpart:
coxph.Lexis <- 
function( Lx, # Lexis object
     formula, # timescale ~ model
        from = preceding(Lx,to), # Exposure ('from' states)
          to = absorbing(Lx)   , # Events ('to' states)
      paired = FALSE,
     verbose = TRUE,
          ... )
{
# Lexis object ?
if( !inherits(Lx,"Lexis") ) stop( "The first argument must be a Lexis object.\n")

# sensible defaults if only one of to and from is missing
if(  missing(from) & !missing(to) ) from <- preceding (Lx,to  )
if( !missing(from) &  missing(to) ) to   <- succeeding(Lx,from)

# name of the dataset
nameLx <- deparse(substitute(Lx))

# work out which transitions are modeled
# first a small utility for annotation
trt <- function( f, t ) paste( f, "->", t, sep="" )
if( paired )
  {
if( length(from) != length(to) ) stop("If 'paired' is TRUE, from and to must have same length!\n")
if( any(from==to) ) stop("If 'paired' is TRUE, entries in from and to must be pairwise different\n")
trnam <- trt( from, to )
  } else {  
tm <- tmat( Lx )[from,to,drop=FALSE]
trnam <- outer( rownames(tm), colnames(tm), trt )[tm>0]
trnam <- trnam[!is.na(trnam)]
  }
# just for formatting explanatory text
onetr <- length( trnam )==1
trprn <- paste( trnam, collapse=", " )
    
# warn if a potentially silly model is defined
if( any( ts<-table(sapply( strsplit(trprn,"->"), function(x) x[1] ))>1 ) ) warning(
 "NOTE:\nMultiple transitions *from* state '",names(ts[ts>1]),"' - are you sure?",
 "\nThe analysis requested is effectively merging outcome states.", 
 "\nYou may want analyses using a *stacked* dataset - see ?stack.Lexis\n" )

# Correct formula?
if( length(formula) != 3 )
    stop("'formula' must be a 2-sided formula, with the l.h.s. the timescale")

# Is the l.h.s. a timescale?
ts <- as.character( formula[2] )
if( !(ts %in% (tms<-timeScales(Lx))) ) 
  stop( "l.h.s. of formula must be a timescale; one of:\n", tms, "\n" )

# What are the 'from' states
from <- levels( factor(Lx$lex.Cst) )
    
# construct a Surv response object, and note that we want the possibility
# of transitions to transient states, hence the lex.Xst != lex.Cst 
Sobj <- Surv( Lx[,ts], 
              Lx[,ts]+Lx$lex.dur,
              trt( Lx$lex.Cst, Lx$lex.Xst ) %in% trnam )
#              Lx$lex.Xst %in% to &
#              Lx$lex.Xst != Lx$lex.Cst )

# Tell what we intend to and then do it    
if( verbose ){
cat( "survival::coxph analysis of Lexis object ", nameLx,
     ":\nRates for", if(  onetr ) " the", " transition",
                     if( !onetr ) "s", " ", trprn,
     "\nBaseline timescale: ", ts, sep="" )
             }

mod <- coxph( as.formula( paste( "Sobj", 
                                 as.character(formula[3]),
                                 sep="~") ), 
              data = Lx, ... )

# Add an explanatory attribute
attr( mod, "Lexis" ) <- list( data=nameLx,
                             trans=trnam,
                           formula=formula )     
class( mod ) <- c( "coxph.lex", class(mod) )    
mod
}
