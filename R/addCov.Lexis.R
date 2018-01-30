# The addCov method
addCov <- function (x, ...) UseMethod("addCov")

addCov.default <-
addCov.Lexis <-
function( Lx,
        clin,
   timescale = 1,
       exnam,
         tfc = "tfc",
   addScales = FALSE )
{
# Function to add clinically measured covariates to a Lexis object
    
# The point is to cut the Lexis diagram at the examination dates
# and subsequently add the clinical records
# ...but first the usual cheking of paraphernalia

if( !inherits(Lx  ,"Lexis") ) stop( "Lx must be a Lexis object.\n" )
if(  inherits(clin,"Lexis") ) stop( "clin cannot be a Lexis object.\n" )
    
# Is the timescale argument a timescale in Lx and is it a variable in clin?    
ts <- if( is.numeric(timescale) ) timeScales( Lx )[timescale] else timescale
if( !( ts %in% timeScales(Lx) ) )
    stop( "timescale argument (", ts, ") must be one of timescales in in the Lexis object ",
          deparse(substitute(Lx)),":", timeScales(Lx), ".\n" )

clin.nam <- deparse(substitute(clin))
if( !( ts %in% names(clin) & "lex.id" %in% names(clin) ) )
    stop( "'lex.id' and timescale '", ts, "' must be a variables in the clin object ",
          clin.nam, "\n" )

# variables to merge by
mvar <- c( "lex.id", ts )

# order clin to get the possible construction of examination names ok
clin <- clin[order(clin$lex.id,clin[,ts]),]

# check that examination dates are unique within persons
if( any( dd <- duplicated(clin[,c("lex.id",ts)]) ) )
    {
  warning( "Examination dates must be unique within persons\n",
           sum(dd), " records with duplicate times from clin object ", clin.nam, 
           " excluded.")
  clin <- clin[!dd,]
    }
    
# the variable holding the name of the examination
if( missing(exnam) ) exnam <- "exnam"
# and if it is not there, construct it
if( !(exnam %in% names(clin)) )
    clin[,exnam] <- paste( "ex",
                           ave( clin$lex.id,
                                clin$lex.id,
                                FUN = function(x) cumsum(x/x) ),
                           sep="" )

# Add copy of the time of examination to be carried forward
clin[,tfc] <- clin[,ts]
    
# Clinical variables to be merged in --- note we take examination date
# and name as a cinical variable too 
cvar <- setdiff( names(clin), mvar )

# A data frame of cutting times
cfr <- data.frame( lex.id = clin$lex.id,
                      cut = clin[,ts],
                new.state = clin[,exnam] )

# Now cut Lx --- this is really inefficient
mc <- Lx
for( st in levels(cfr$new.state) )
mc <- cutLexis( mc,
               cut = cfr[cfr$new.state==st,],
         timescale = ts,
  precursor.states = NULL,
         new.scale = addScales )
    
# Merge in states from the original object mx, but take attributes from mc
mx <- Lx[,mvar]
mx$org.Cst <- Lx$lex.Cst
mx$org.Xst <- Lx$lex.Xst
mx <- merge( mc, mx, by = mvar, all.x = TRUE, sort = TRUE )    
    
# Complete the state variables    
( wh <- which(is.na(mx$org.Cst)) )
mx$org.Cst[wh] <- na.locf( mx$org.Cst, nx.rm=FALSE )[wh]
mx$org.Xst[wh] <- na.locf( mx$org.Xst, nx.rm=FALSE )[wh]
# but - oops - the earlier lex.Xst should be as lex.Cst
mx$org.Xst[wh-1] <- mx$org.Cst[wh-1]
# overwrite the useless ones and get rid of the intermediate
mx$lex.Cst <- mx$org.Cst
mx$lex.Xst <- mx$org.Xst
wh.rm <- match( c("org.Cst","org.Xst"), names(mx) )
mx <- mx[,-wh.rm]
    
# Merge in the clinical variables
mx <- merge( mx, clin, by=mvar, all.x=TRUE, sort=TRUE )

# And carry them forward within each lex.id
    
# locf within each person (should be easier in data.table)
locf.df <- function( df ) as.data.frame( lapply( df, na.locf, na.rm=FALSE ) )

# ave does not like character variables so we convert to factors
wh <- which( sapply( mx[,cvar], is.character ) )
for( j in wh ) mx[,cvar[j]] <- factor( mx[,cvar[j]] )
# then we can carry forward
mx[,cvar] <- ave( mx[,cvar], mx$lex.id, FUN=locf.df )

# Finally update the time from last clinical measurement
mx[,tfc] <- mx[,ts] - mx[,tfc]

# Add as a time-scale
if( addScales )
  {  
new.scales <- setdiff( timeScales(mx), timeScales(Lx) )
op <- options(warn = (-1)) # suppress warnings
mx[,tfc] <- apply( mx[,new.scales,drop=FALSE], 1, min, na.rm=TRUE )    
options(op) # reset the previous value
attr( mx, "time.scales") <- c( attr( mx, "time.scales"), tfc ) 
attr( mx, "time.since" ) <- c( attr( mx, "time.since"), "" )
brt <- list( x=NULL ) ; names( brt ) <- tfc
attr( mx, "breaks") <- c( attr( mx, "breaks"), brt ) 
  }
    
# Done! - well order first
mx[order(mx[,"lex.id"],mx[,timeScales(mx)[1]]),]
}
