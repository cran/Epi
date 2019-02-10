# Acronyms:
# dop: date of purchase
# dpt: dose per time
# amt: amount
# dur: duration
# dof: date of follow-up

######################################################################
# A function to compute cumulative doses etc. at the purchase dates
# *and* derive intervals of drug coverage resp no coverage from the
# purchase date, amounts and dosage per time
use.amt.dpt <-
function( purchase, # data frame with dop, amt and dpt
          push.max,
          rm.dose )
{
do.call( "rbind",
lapply( split( purchase, purchase$id ),
        function(set)
        {
        np <- nrow(set)
        if( np==1 ) # return data frame of 2 rows
          {
          dfR <- data.frame( id = set$id,
                            dof = set$dop + c(0,set$amt/set$dpt),
                            amt = c(set$amt,0),
                            dpt = c(set$dpt,0),
                            off = c(FALSE,TRUE),
                        cum.amt = c(0,set$amt),
                        cum.tim = c(0,set$amt/set$dpt) )
          return( dfR )
          }
        set <- set[order(set$dop),]
        # Compute length of exposure periods
        drug.dur  <- set$amt / set$dpt
        # Put the exposed period head to foot
        new.start <- min( set$dop ) + c(0,cumsum(drug.dur[-np]))
        # Move them out so that the start of a period is never earlier than
        # the dop
        exp.start <- new.start + cummax( pmax(set$dop-new.start,0) )
        # Compute the pushes
        push.one  <- exp.start - set$dop
        # Revise them to the maximally acceptable
        push.adj  <- pmin( push.one, push.max )
        # Revise the starting dates of exposure
        exp.start <- exp.start - push.one + push.adj
        # Revise the durations to be at most equal to differences between the
        # revised starting dates
        drug.dur  <- pmin( drug.dur, c(diff(exp.start),Inf) )
        # Compute the end of the intervals
        exp.end   <- exp.start + drug.dur
        # Find out which exposure intervals that are followed by a gap:
        followed.by.gap <- c( exp.start[-1]-exp.end[-np] > 0, TRUE )
        # To facilitate further calculations we add records for the gaps
        # and also a record for the date of drug expiration
        dfR <- rbind( data.frame( id = set$id[1],
                                 dof = exp.start,
                                 amt = set$amt,
                                 dpt = set$dpt ),
                      data.frame( id = set$id[1],
                                 dof = exp.end[followed.by.gap],
                                 amt = 0,
                                 dpt = 0 ) )
        dfR <- dfR[order(dfR$dof),]
        # Logical indicator of intervals off drug
        dfR$off <- (dfR$dpt==0) | (dfR$dof < min(set$dop))
        # Finally compute the cumulative dose and time on drug at the end of
        # the interval using interval length and dpt:
        if( rm.dose ) {
        dfR$cum.amt <- with( dfR, cumsum( c(0, diff(dof)* dpt[-length(dpt)]) ) )
        } else dfR$cum.amt <- cumsum( c( 0, dfR$amt[-nrow(dfR)] ) )
        dfR$cum.tim <- with( dfR, cumsum( c(0, diff(dof)*!off[-length(dpt)]) ) )
        return( dfR )
        } ) )
}

######################################################################
# Compute the cumulative dose at all purcase dates and at the last
# (unknown) future expiry date, computed based on previous
# consumption.  The resulting data frame has one more line per person
# than no. of purchases.
use.only.amt <-
function( purchase,
          pred.win )
{
do.call( "rbind",
lapply( split( purchase, purchase$id ),
        function(set)
        {
        np <- nrow(set)
        if( np==1 ) return( NULL )
        set <- set[order(set$dop),]
        # The points to include in the calculation:
        # All dates after pred.win before last purchase,
        # but at least the last two purchase dates,
        wp <- ( set$dop > pmin( max(set$dop)-pred.win,
                               sort(set$dop,decreasing=TRUE)[2] ) )
        # Cumulative amount consumed at each dop since first
        cum.amt <- cumsum(c(0,set$amt))
        # Average slope to use to project the duration of last purchase
        avg.slp <- diff(range(cum.amt[c(wp,FALSE)]))/
                   diff(range(set$dop[wp]))
        # Purchase dates and the date of last consumption
        dof <- c( set$dop, set$dop[np]+set$amt[np]/avg.slp )
        # Cumulative time since first
        cum.tim <- dof - set$dop[1]
        # The only time without drug is after last dispense expired
        off <- rep( c(FALSE,TRUE), c(length(dof)-1,1) )
        return( data.frame( id = set$id[1],
                           dof = dof,
                       cum.amt = cum.amt,
                       cum.tim = cum.tim,
                           off = off ) ) # will be used in the subsequent code
        } ) )
}

######################################################################
# Function to interpolate the drug exposures to the dates of FU.
dex.int <-
function( set, breaks, lags, lag.dec )
{
  dof <- NULL
  # If only one record return null
  if( nrow(set) < 2 ) return( NULL )
  # All values of these are identical within each set (=person)
  doe <- set$doe[1]
  dox <- set$dox[1]
  # The first date of drug exposure according to the assumption
  doi <- min(set$dof)
  # Breakpoints and the entry end exit dates are the dates we care for
  # - but only within breaks
  doin <- max( min(breaks), doe )
  doex <- min( max(breaks), dox )
  breaks <- sort( unique( c(breaks,doin,doex) ) )
  # only the breaks inside the follow-up interval
  xval <- breaks[breaks>=doin & breaks<=doex]
  if( length(xval) == 0 ) return( NULL )
  # Merge a) the purchase dates (set$dof) and
  #    indicator being off drug (set$off) and
  # the date a person goes off drug, doff (set to NA for remaining records)
  # with b) the break dates (which is where we want things computed)
  # Note we merge on the variable dof and have the data frame with
  # dof=xval as the >first< so that values of xval will be in the
  # resulting dfr$dof even if dof-values in set are almost equal til
  # an xval value. Ensures that sum(dfr$dof %in% xval)==length(xval)   
  dfr <-  merge( data.frame( id = set$id[1],
                            dof = xval ),
                 data.frame( set[,c("id","dof","off")],
                             doff = ifelse(set$off,set$dof,NA) ),
                 all=TRUE )
  # carry the off drug indicator forward to all break dates
  dfr$off <- zoo::na.locf( dfr$off, na.rm=FALSE )
  # possible NAs are before any drug purchase hence off-drug
  dfr$off <- ifelse( is.na(dfr$off), TRUE, dfr$off )
  # carry date of going off drug forward, but only for the off period
  dfr$doff <- zoo::na.locf( dfr$doff, na.rm=FALSE ) * ifelse(dfr$off,1,NA)
  # time from cessation can now be computed
  dfr$tfc <- dfr$dof - dfr$doff
  # time from initiation of drug
  dfr$tfi  <- pmax( 0, dfr$dof-doi )
  # restrict to the desired timepoints
  dfr <- subset( dfr, dof %in% xval )
  dfr <- dfr[!duplicated(dfr$dof),]
  # linear interpolation of the cumulative dose and time from the
  # purchase data (set)
  dfr$cdos <- approx( set$dof, set$cum.amt, xout=xval, rule=2 )$y
  dfr$ctim <- approx( set$dof, set$cum.tim, xout=xval, rule=2 )$y
  # the same for the desired lags
  for( lg in lags )
     dfr[,paste( "lag.pre",
                 formatC(lg,format="f",digits=lag.dec),
                 sep="" )] <-
     approx( set$dof, set$cum.amt, xout=xval-lg, rule=2 )$y
  dfr$doff <- zoo::na.locf( dfr$doff, na.rm=FALSE )
  dfr$tfc[is.na(dfr$tfc)] <- 0
  dfr$dur <- c( diff(dfr$dof), NA )
  return( dfr[-nrow(dfr),] )
  }

######################################################################
# The function that ties it all together
gen.exp <-
function( purchase,  id="id" , dop="dop", amt="amt", dpt="dpt",
                fu, doe="doe", dox="dox",
            breaks,
           use.dpt = ( dpt %in% names(purchase) ),
          push.max = Inf,
           rm.dose = FALSE,
              lags = NULL,
           lag.dec = 1,
           lag.pre = "lag.",
          pred.win = Inf )
{
# to aviod a NOTE from R CMD check
dof <- NULL
# Make sure that the data frames have the right column names
wh <- match( c(id,dop,amt), names(purchase) )
if( any( is.na(wh) ) ) stop("Wrong column names for the purchase data frame")
names( purchase )[wh] <- c("id","dop","amt")

# Allow dpt to be entered as numerical scalar common for all records
if( use.dpt )
  if( is.numeric(dpt) )
    {
    if( length(dpt) > 1 ) stop("If dpt is numeric it must have length 1\n",
        "otherwise it must be the name of the dpt variable in the dataset")
    purchase$dpt <- dpt
    }
  else names( purchase )[match(dpt,names(purchase))] <- "dpt"

wh <- match( c(id,doe,dox), names(fu) )
if( any( is.na(wh) ) ) stop("Wrong column names for the follow-up data frame")
names( fu )[wh] <- c("id","doe","dox")

if( use.dpt ) { tmp.dfr <- use.amt.dpt( purchase,
                                        push.max = push.max,
                                        rm.dose = rm.dose )
       } else { tmp.dfr <- use.only.amt( purchase,
                                         pred.win = pred.win ) }

# Having done what needs to be done about the purchase records we turn
# to the follow-up records, but we first need to accommodate the fact
# that there might be more fu records per person:

# Generate record no within id in fu (ordering in fu is immaterial)
fu$no <- ave( fu$id, fu$id, FUN=function(x) 1:length(x) )
# Set up the object to collect the resulting follow-up
res <- NULL

# loop through subsets of fu with different ids
for( ni in 1:max(fu$no) )
   {
# Merge the follow-up period for the persons to the correponding
# exposure records (hence the all.y=T)
tm.dfr <- merge( tmp.dfr, fu[fu$no==ni,], by="id", all.y=TRUE )

# Interpolate to find the cumulative doses at the dates in the vector breaks
res.dfr <-
do.call( "rbind",
         lapply( split( tm.dfr, tm.dfr$id ),
                 dex.int,
                 breaks, lags, lag.dec ) ) # end of do.call
res <- rbind( res, res.dfr ) # append to result from other fu$no
   } # end of for( ni in 1:max(fu$no) )
    
var.order <- c("id","dof","dur","off","doff","tfc","tfi","ctim","cdos",
               names(res)[grep("lag.prefix",names(res))] )
res <- res[order(res$id,res$dof),var.order]
names(res) <- gsub( "lag.prefix", lag.pre, names(res) )
res
}
# end of gen.ex
