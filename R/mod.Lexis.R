modLexis <- 
function( Lx, resp, formula, xpos, link="log", scale, verbose=TRUE, ..., model )
{
# A common wrapper for glm and gam modeling of Lexis FU
  
# is this a Lexis object ?
if( !inherits(Lx,"Lexis") ) stop( "The first argument must be a Lexis object.\n")
nameLx <- deparse(substitute(Lx))
# Beginning of a new feature with a countmultiplier of the transitions considered 
if( !("lex.N" %in% names(Lx)) ) Lx$lex.N <- 1
    
# check that events are actual levels of lex.Xst 
if( is.numeric(resp) ) resp <- levels( Lx$lex.Xst )[resp]
wh <- match( resp, levels(Lx$lex.Xst) )
if( any(is.na(wh)) ) stop("'resp' must be a subset of: '", 
                     paste(levels(Lx$lex.Xst), collapse="','", sep=""), "'\n" )  

# is xpos supplied?
if( missing(xpos) ) {
  xpos <- levels( factor(Lx$lex.Cst) )
  } else {
# check that xpos are actual levels of lex.Cst 
if( is.numeric(xpos) ) xpos <- levels( Lx$lex.Cst )[xpos]
wh <- match( xpos, levels(Lx$lex.Cst) )
if( any(is.na(wh)) ) stop("'xpos' must be a subset of: '", 
                     paste(levels(Lx$lex.Cst), collapse="','", sep=""), "'\n" )  
Lx <- Lx[Lx$lex.Cst %in% xpos,]                     
  } 
                     
# construct the model formula - note that we want the possibility of
# transitions to transient states, hence the lex.Xst != lex.Cst
if( length(formula) != 2 ) stop( "formula must be a one-sided formula")
form <- cbind( (Lx$lex.Xst %in% resp &
                Lx$lex.Xst != Lx$lex.Cst)*Lx$lex.N,
               Lx$lex.dur ) ~ 1
form[3] <- formula[2]
xpos <- levels( factor(Lx$lex.Cst) ) # only levels present in lex.Cst

# Scaling
Lx$lex.dur <- Lx$lex.dur/scale
    
# Tell what we intend to and then do it
if( verbose ){
cat( deparse(substitute(model)),
     " Poisson analysis of Lexis object ", nameLx, " with ", link, " link",
     ":\n Transition rates from '", paste( xpos, collapse="','"), 
                          "' to '", paste( resp, collapse="','"), "'",
     if( scale!=1 ) paste(" scaled by", scale ), "\n", sep="" )
             }
    
# Fit the model
mod <- model( form, family = poisreg(link=link), data = Lx, ... )

# An explanatory attribute
attr( mod, "Lexis" ) <- list( Exposure=xpos, Events=resp, scale=scale )     
mod
}

# Here are the actual functions of interest

glm.Lexis <- 
function( Lx, resp, formula, xpos, link="log", scale=1    , verbose=TRUE   , ... ) {
modLexis( Lx, resp, formula, xpos, link=link , scale=scale, verbose=verbose, ..., model=stats::glm ) }

gam.Lexis <- 
function( Lx, resp, formula, xpos, link="log", scale=1    , verbose=TRUE   , ... ) {
modLexis( Lx, resp, formula, xpos, link=link , scale=scale, verbose=verbose, ..., model= mgcv::gam ) }

# and here the coxph counterpart:

coxph.Lexis <- 
function( Lx, # Lexis object 
        resp, # Events ('to' states)
     formula, # timescale ~ model
        xpos, # exposure states ('from' states)
     verbose = TRUE,
          ... )
{
# Lexis object ?
if( !inherits(Lx,"Lexis") ) stop( "The first argument must be a Lexis object.\n")
nameLx <- deparse(substitute(Lx))

# check levels    
if( is.numeric(resp) ) resp <- levels( Lx$lex.Xst )[resp]
wh <- match(resp,levels(Lx$lex.Xst))
if( any(is.na(wh)) ) stop("'resp' must be a subset of: '", 
                     paste(levels(Lx$lex.Xst), collapse="','", sep=""), "'\n" )  

# is xpos supplied?
if( !missing(xpos) )
  {
# check that xpos are actual levels of lex.Cst 
if( is.numeric(xpos) ) xpos <- levels( Lx$lex.Cst )[xpos]
wh <- match( xpos, levels(Lx$lex.Cst) )
if( any(is.na(wh)) ) stop("'xpos' must be a subset of: '", 
                     paste(levels(Lx$lex.Cst), collapse="', '", sep=""), "'\n" )  
Lx <- Lx[Lx$lex.Cst %in% xpos,]                     
  } 

# Correct formula?
if( length(formula) != 3 ) stop("'formula' must be a 2-sided formula, with the l.h.s. the timescale")

# Is the l.h.s. a timescale?
ts <- as.character( formula[2] )
if( !(ts %in% (tms<-timeScales(Lx))) ) 
  stop( "l.h.s. of formula must be a timescale; one of:\n", tms, "\n" )

# What are the 'from' states
xpos <- levels( factor(Lx$lex.Cst) )
    
# construct a Surv response object, and note that we want the possibility
# of transitions to transient states, hence the lex.Xst != lex.Cst 
Sobj <- Surv( Lx[,ts], 
              Lx[,ts]+Lx$lex.dur,
              Lx$lex.Xst %in% resp &
              Lx$lex.Xst != Lx$lex.Cst )

# Tell what we intend to and then do it    
cat( "survival::coxph analysis of Lexis object ", nameLx, " using timescale ", ts,
     ":\nTransition rates from '", paste( xpos, collapse="','"), 
                         "' to '", paste( resp, collapse="','"),  "'\n", sep="" )
mod <- coxph( as.formula( paste( "Sobj", 
                                 as.character(formula[3]),
                                 sep="~") ), 
              data = Lx, ... )

# An explanatory attribute
attr( mod, "Lexis" ) <- list( Exposure=xpos, Events=resp, Timescale=ts )
mod
}
