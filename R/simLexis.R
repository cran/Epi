# First the utility functions

cummid <-
function( x, time.pts=1:length(x) )
{
# Computes the cumulative area under a curve vith values x at time.pts
cumsum( c(0, (x[-1]-diff(x)/2)*diff(time.pts) ) )
}

sim1 <-
function( rt, init, time.pts )
{
# Simulates a single transition time and state based on the dataframe
# rt with columns lex.id and timescales. Each row in rt is the id,
# followed by the set of estimated transition rates to the different
# states reachable from the current one.
ci <- apply( rt[,-1,drop=FALSE], 2, cummid, time.pts )
tt <- uu <- -log( runif(ncol(ci)) )
for( i in 1:ncol(ci) ) tt[i] <- approx( ci[,i], time.pts, uu[i], rule=2 )$y
# Note this resulting data frame has 1 row
data.frame( lex.id  = rt[1,1],
            lex.dur = min(tt,na.rm=TRUE),
            lex.Xst = factor( if( min(tt)<max(time.pts) ) colnames(ci)[tt==min(tt)]
                              else NA, levels=levels(init$lex.Cst) ) )
}

simX <-
function( nd, init, Tr, time.pts )
{
# Simulation is done from the data frame nd, in chunks of starting
# state, lex.Cst. This is necessary because different states have
# different (sets of) exit rates. Therefore, this function simulates
# for a set of persons from the same starting state.
np <- length( time.pts )
nr <- nrow( nd )
if( nr==0 ) return( NULL )

# The as.character below is necessary because indexing by a factor
# by default is by the number of the level, and we are not indexing by
# this, but by components of Tr which just happens to have  names that
# are a subset of the levels of lex.Cst.
cst <- as.character( unique(nd$lex.Cst) )
if( length(cst)>1 ) stop( "More than one lex.Cst present.\n" )

# Expand each person by the timepoints
nx <- nd[rep(1:nrow(nd),each=np),]
nx[,timeScales(init)] <- nx[,timeScales(init)] + rep(time.pts,nr)
nx$lex.dur <- 1

# Make a dataframe with predicted rates for each of the transitions
# out of this state for these times
rt <- data.frame( lex.id=nx$lex.id )
for( i in 1:length(Tr[[cst]]) ) rt <- cbind( rt, exp(predict(Tr[[cst]][[i]],newdata=nx)) )
names( rt )[-1] <- names( Tr[[cst]] )

# Then find the transition time and exit state for each person:
xx <- match( c("lex.dur","lex.Xst"), names(nd) )
if( any( !is.na(xx) ) ) nd <- nd[,-xx[!is.na(xx)]]
merge( nd, do.call( "rbind", lapply( split(rt,rt$lex.id), sim1, init, time.pts ) ), by="lex.id" )
}

get.next <-
function( sf, init, tr.st )
{
# Procduces an initial Lexis object for the next simulation for those
# who have ended up in a transient state.
# Note that this exploits the existance of the "time.since" attribute
# for Lexis objects and assumes that a character vector naming the
# transient states is supplied as argument.
if( nrow(sf)==0 ) return( sf )
nxt <- sf[sf$lex.Xst %in% tr.st,]
if( nrow(nxt) == 0 ) return( nxt )
nxt[,timeScales(init)] <- nxt[,timeScales(init)] + nxt$lex.dur
wh <- attr( init,"time.since" )
for( i in 1:length(wh) )
   if( wh[i] != "" ) nxt[nxt$lex.Xst==wh[i],timeScales(init)[i]] <- 0
nxt$lex.Cst <- nxt$lex.Xst
return( nxt )
}

chop.lex <-
function( obj, cens )
{
# A function that chops off all follow-up beyond cens since entry for
# each individual
zz <- entry( obj, 1, by.id=TRUE )
ww <- merge( obj, data.frame( lex.id=as.numeric(names(zz)),
                                cens=zz+cens ) )
ww <- ww[ww[,timeScales(obj)[1]] < ww$cens,]
x.dur <- pmin( ww$lex.dur, ww[,"cens"]-ww[,timeScales(obj)[1]] )
ww$lex.Xst[x.dur<ww$lex.dur] <- ww$lex.Cst[x.dur<ww$lex.dur]
ww$lex.dur <- pmin( x.dur, ww$lex.dur )
ww
}

simLexis <-
function( Tr, # List of lists of glm objects
        init, # Lexis objects of persons to simulate. Must have the
              # same attributes as the original object, in particular
              # "time.scales" and "time.since".
    time.pts = 0:50/2, # Points where rates are computed in the
                       # simulation
           N = 1, # How many persons should each line in
                  # init represent?
      lex.id = 1:(N*nrow(init)), # What should be the ids of the simulated persons?
        type = "glm-mult" # Not used currently
         )
{
# Expand the input data frame using N and put in lex.id
foo <- lex.id
init <- init[rep(1:nrow(init),each=N),]
init$lex.id <- foo

# Fix attributes
if( is.null( nts <- attr(init,"time.scales") ) )
  stop( "No time.scales attribute for init" )
if( is.null( attr(init,"time.since") ) )
  {
  attr(init,"time.since") <- rep( "", nts )
  cat( "WARNING:\n
       'time.since' attribute set, which means that you assume that\n
        none of the time scale represent time entry to a state." )
  }
# Convenience constants
np <- length( time.pts )
tr.st <- names( Tr )

# The first set of sojourn times in the initial states
sf <- do.call( "rbind", lapply( split(init,init$lex.Cst), simX, init, Tr, time.pts ) )

# Then we must update those who have ended in transient states
# and keep on doing that till all are in absorbing states or censored
nxt <- get.next( sf, init, tr.st )
while( nrow(nxt) > 0 )
{
nx <- do.call( "rbind", lapply( split(nxt,nxt$lex.Cst), simX, init, Tr, time.pts ) )
sf <- rbind( sf, nx )
nxt <- get.next( nx, init, tr.st )
}

# Doctor lex.Xst for the censored, and supply attributes
sf$lex.Xst[is.na(sf$lex.Xst)] <- sf$lex.Cst[is.na(sf$lex.Xst)]

# Finally, nicely order the output by persons, then times and states
nord <- match( c( "lex.id", timeScales(sf),
                  "lex.dur",
                  "lex.Cst",
                  "lex.Xst" ), names(sf) )
noth <- setdiff( 1:ncol(sf), nord )
sf <- sf[order(sf$lex.id,sf[,timeScales(init)[1]]),c(nord,noth)]
rownames(sf) <- NULL
attr( sf, "time.scales" ) <- attr( init, "time.scales" )
attr( sf, "time.since"  ) <- attr( init, "time.since" )
chop.lex( sf, max(time.pts) )
}

nState <-
function ( obj,
            at,
          from,
     time.scale = 1 )
{
# counte the number of persons in each state of the Lexis object 'obj'
# at the times 'at' from the time 'from' in the time scale
# 'time.scale'

# Determin timescales and absorbing and transient states
tmsc <- Epi:::check.time.scale(obj,time.scale)
TT <- tmat(obj)
absorb <- rownames(TT)[apply(!is.na(TT),1,sum)==0]
transient <- setdiff( rownames(TT), absorb )

# Expand each record length(at) times
tab.frm <-
    obj[rep(1:nrow(obj),each=length(at)),c(tmsc,"lex.dur","lex.Cst","lex.Xst")]

# Stick in the correponding times on the chosen time scale
tab.frm$when <- rep( at, nrow(obj) ) + from

# For transient states keep records that includes these points in time
tab.tr <- tab.frm[tab.frm[,tmsc]                 <= tab.frm$when &
                  tab.frm[,tmsc]+tab.frm$lex.dur >  tab.frm$when,]
tab.tr$State <- tab.tr$lex.Cst

# For absorbing states keep records where follow-up ended before
tab.ab <- tab.frm[tab.frm[,tmsc]+tab.frm$lex.dur <= tab.frm$when &
                  tab.frm$lex.Xst %in% absorb,]
tab.ab$State <- tab.ab$lex.Xst

# Make a table using the combination of those in transient and
# absorbing states.
with( rbind( tab.ab, tab.tr ), table( when, State ) )
}

pState <-
function( nSt, perm=1:ncol(nSt) )
{
# Compute cumulative proportions of persons across states in order
# designate by 'perm'
tt <- t( apply( nSt[,perm], 1, cumsum ) )
tt <- sweep( tt, 1, tt[,ncol(tt)], "/" )
class( tt ) <- c("pState","matrix")
tt
}

plot.pState <-
function( x,
        col = rainbow(ncol(x)),
     border = "transparent",
       xlab = "Time",
       ylab = "Probability", ... )
{
# Function to plot cumulative probabilities along the time scale.

# Just for coding convenience when plotting polygons
pSt <- cbind( 0, x )
matplot( as.numeric(rownames(pSt)), pSt, type="n",
         ylim=c(0,1), yaxs="i", xaxs="i",
         xlab=xlab, ylab=ylab, ... )
for( i in 2:ncol(pSt) )
   {
   polygon( c(    as.numeric(rownames(pSt)) ,
              rev(as.numeric(rownames(pSt))) ),
            c(    pSt[,i  ],
              rev(pSt[,i-1]) ),
            col=col[i-1], border=border[i-1], ... )
   }
}
