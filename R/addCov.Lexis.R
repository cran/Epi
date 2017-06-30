# The factorize method
addCov <- function (x, ...) UseMethod("addCov")

addCov.default <-
addCov.Lexis <-
function( Lx, clin,
          timeScale = 1,
          addScales = FALSE,
          exnam,
          tfc = "tfc" )
{
# Function to merge clinical information into a multistate object
if( !inherits(Lx  ,"Lexis") ) stop( "Lx must be a Lexis object.\n" )
if(  inherits(clin,"Lexis") ) stop( "clin cannot be a Lexis object.\n" )
    
# Is the timeScale argument a timescale and is it in clin?    
ts <- if( is.numeric(timeScale) ) timeScales( Lx )[timeScale] else timeScale
if( !( ts %in% timeScales(Lx) ) )
    stop( "timeScale (", ts, ") must be in the Lexis object ",
          deparse(substitute(Lx)), "\n" )
if( !( ts %in% names(clin) ) )
    stop( "timeScale (", ts, ") must be a variable in clin object ",
          deparse(substitute(clin)), "\n" )

# variables to merge by
mvar <- c( "lex.id", ts )

# order clin to get the possible construction of examination names ok
clin <- clin[order(clin$lex.id,clin[,ts]),]

# the variable holding the name of the examination
if( missing(exnam) ) exnam <- "exnam"
# if it is not there, construct it
if( !(exnam %in% names(clin)) )
    clin[,exnam] <- paste( "ex",
                           ave( clin$lex.id,
                                clin$lex.id,
                                FUN = function(x) cumsum(x/x) ),
                           sep="." )

# check examination names are unique within persons
if( any( table(clin$lex.id,clin[,exnam])>1 ) )
  stop( "Examination names must be unique within persons\n" )

# Add copy of the time of examination to be carried forward
clin[,tfc] <- clin[,ts]
    
# Clinical variables to be merged in --- note we take examination date in
cl.vars <- setdiff( names(clin), mvar )

# Merge in the clinical variables at the correct times
# In most cases this will just result in a frame with
# nrow(Lx)+nrow(clin) rows
mx <- merge( Lx, clin, by=mvar, all=T )

# order rows by time within each person 
mx <- mx[order(mx$lex.id,mx[,ts]),]

# locf within each person (should be easier in data.table)
locf.df <- function( df ) as.data.frame( lapply( df, na.locf, na.rm=FALSE ) )
# ave does not like character variables so we convert to character
wh <- which( sapply( mx[,cl.vars], is.character ) )
for( j in wh ) mx[,cl.vars[j]] <- factor( mx[,cl.vars[j]] )
mx[,cl.vars] <- ave( mx[,cl.vars], mx$lex.id, FUN=locf.df )

# a copy of the original state variables (essentially renaming them)
mx$org.Cst <- mx$lex.Cst
mx$org.Xst <- mx$lex.Xst

# Carry state info forward --- this is necessarily within a person
# --- we later remove records before entry for each person, so no worries
( wh <- which(is.na(mx$org.Cst)) )
mx$org.Cst[wh] <- na.locf( mx$org.Cst, nx.rm=FALSE )[wh]
mx$org.Xst[wh] <- na.locf( mx$org.Xst, nx.rm=FALSE )[wh]
# but - oops - the earlier lex.Xst should be as lex.Cst
mx$org.Xst[wh-1] <- mx$org.Cst[wh-1]
    
# need these as charcter variables not factors
mx$org.Cst <- as.character(mx$org.Cst)
mx$org.Xst <- as.character(mx$org.Xst)

# Cut at the examnation dates to get the timescales and lex.dur correct
cfr <- data.frame( lex.id = clin$lex.id,
                      cut = clin[,ts],
                new.state = clin[,exnam] )
mc <- Lx
# This is where we assume that examinatin names are unique in persons
for( st in unique(as.character(cfr$new.state)) )
mc <- cutLexis( mc, cut = cfr[cfr$new.state==st,],
              timescale = ts,
       precursor.states = levels(Lx) )
mc$lex.Cst <- as.character(mc$lex.Cst)
mc$lex.Xst <- as.character(mc$lex.Xst)

# meaningless to merge two Lexis frames; mc is the one carrying the
# Lexis struture of interest, so mx is the slave (hence org.Cst, org.Xst):
class( mx ) <- "data.frame"

Lr <- merge( mc, mx[,c(mvar,cl.vars,"org.Cst","org.Xst")], by=mvar, all.x=TRUE )
Lr$lex.Cst <- Lr$org.Cst
Lr$lex.Xst <- Lr$org.Xst

# Time since last covariate time
Lr[,tfc] <- Lr[,ts] - Lr[,tfc]

# Make lex.Cst and lex.Xst to proper factors again
Lr <- Relevel( Lr[,-match(c("org.Cst","org.Xst"),names(Lr))] )

# Useful function --- later 
order.Lexis <- function( x, ... ) order( x$lex.id, x[,timeScales(x)[1]], ... )
 sort.Lexis <- function( x, ... ) x[order.Lexis( x, ... ),]    
# Order the data.frame
sort.Lexis( Lr )
}
