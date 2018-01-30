mcutLexis <-
function( L0, # A Lexis object
       timescale = 1,    # the time scale referred to by L0[,wh]
              wh,        # indices/names of columns holding dates of state entries (events)
      new.states = NULL, # Names of the event types (states)
precursor.states = NULL,
      seq.states = TRUE, # Should state names reflect ordering of events
      new.scales = NULL, # Time-scales referring to time since
    ties.resolve = FALSE # Are tied event times accepted?
        )
{
### we rely on referring to the timescale and event time variables by name
if( is.numeric(timescale) ) timescale <- timeScales(L0)[timescale]    
if( is.numeric(wh) ) wh <- names(L0)[wh]    

### don't be silly
if( length(wh)==1 )
# return( docut( L0, osv ) ) # old cutLexis should be absorbed here
  stop( "mcutLexis not needed for one type of event - use cutLexis\n" )

### states    
if( is.null(new.states) )
  {
  new.states <- wh
  cat( "NOTE: Name of new states set to\n", new.states )
  }
if( length(wh) != length(new.states) ) 
  stop( "wh and new.states must have same length, but lengths are",
        "wh:", length(wh), "and new.states:", length(new.states), "\n" )

### timescales
# either all or none    
if( is.logical(new.scales) )
  if( any( new.scales ) )
  {  
  new.scales <- paste( "tf", new.states, sep="" )
  cat( "NOTE: new.scales set to: ", new.scales, "\n" )
  }
if( is.character(new.scales) & length(new.scales) != length(wh) )
  {
  new.scales <- paste( "tf", new.states, sep="" )
  warning( "new.scales not of same length as wh. Set to: ",
           new.scales, "\n" )
  }
if( is.character(new.scales) & length(intersect(new.states,timeScales(L0))) )
  stop( "Names of new time scales must be different from names of timescales:\n",
        timeScales(L0) )  
    
### Tied transition times untied    
has.ties <- any( wh.tied <- apply( L0[,wh], 1,
                            function(x) any(diff(sort(x[!is.na(x)]))==0) ) )
if( has.ties & is.logical(ties.resolve) & !ties.resolve )
  stop( "Tied event times not allowed with ties.resolve=FALSE:\n",
        "there were", length(wh.tied), "records with tied event times.")
if( has.ties & is.logical(ties.resolve) & ties.resolve ) ties.resolve = 1/100
if( has.ties & is.numeric(ties.resolve) )
  {
  rnd <- L0[wh.tied,wh]*0
  rnd[,] <- runif(rnd,-1,1) * ties.resolve
  L0[wh.tied,wh] <- L0[wh.tied,wh] + rnd
  cat( "NOTE: ", length(wh.tied),
       "records with tied events times resolved.\n",
       "Results only reproducible if the seed for the random number generator is set.")
  }
# End of checks

# The object to return initiated as NULL
Lcut <- NULL

# Utility function returning sequences of ocurrences as paste of numbers
NAorder <- 
function (x) 
    {
    oo <- order(x, na.last = T)
    on <- (1:length(oo))[oo]
    on[is.na(x[oo])] <- NA
    paste(on[!is.na(on)], collapse = "-")
    }
    
# where do the different sequences of events actually occur in data
L0$whseq <- apply( L0[,wh], 1, NAorder )

# Loop through the actually occurring orders of event occurrences
for( sq in unique(L0$whseq) )
{ 
# Persons with none of the events occurring transferred to result
if( sq=="" ) Lcut <- rbind( Lcut, L0[L0$whseq=="",] )
else {
    
# Extract the subset of persons with a given sequence of events
Ltmp <- L0[L0$whseq==sq,]

# The numerical sequence of states (refer to the elements of wh)
ost <- as.numeric( strsplit( sq, "-" )[[1]] )
nxst <- ""
prst <- precursor.states
for( cs in ost )
   {
 nxst <- ifelse( cs==ost[1],  new.states[cs],
                 paste( nxst, new.states[cs], sep="-" ) )
 Ltmp <- cutLexis( Ltmp, cut = Ltmp[,wh[cs]],
                   timescale = timescale,
                   new.state = nxst,
            precursor.states = prst )
 # include the created state among the precursor states for next cut
 prst <- c(prst,nxst)
   } # end of for loop through events in this sequence (cs)
    
# Attach it to the end of the Lexis object
Lcut <- rbind( Lcut, Ltmp )
     } # end of the else clause

} # end of for loop through sequences (sq)

# Do we want the sequences, the unordered set of previous events or
# just the current one:
old.seq <- seq.states
if( is.logical(seq.states) ) seq.states <- ifelse( seq.states, "s", "u" )
if( is.character(seq.states) ) seq.states <- tolower( substr(seq.states,1,1) )
if( !(seq.states %in% c("s","o","u","l","c")) )
    stop( "What do you mean by seq.states=", old.seq,
          "? - it should abbreviate to one of s, o, u, l or c \n")
# Unordered or last (current) states    
if( seq.states %in% c("u","l","c") ) 
  {
  # Each list element is a vector of states visited
  slvl <- strsplit( levels( Lcut ), "-" )
  # merge those that have the same elements or take the last
  rlvl <- if( seq.states=="u" ) { sapply( lapply( slvl, sort ), paste, collapse="+" )
                           } else sapply( slvl, function(x) x[length(x)] )
  # Relevel the states    
  levels( Lcut$lex.Cst ) <-
  levels( Lcut$lex.Xst ) <- rlvl
  }

# Did we ask for timescales as time since events?
if( !is.null(new.scales) )
  {
  # insert columns for the new time scales
  Lcut <- Lcut[,c(rep("whseq",length(new.scales)),names(Lcut))]
  names( Lcut )[1:length(new.scales)] <- new.scales
  for( i in 1:length(wh) )
     Lcut[,i] <- ifelse( Lcut[,timescale] - Lcut[,wh[i]] < 0,
                         NA,
                         Lcut[,timescale] - Lcut[,wh[i]] )
  # set attributes
  attr( Lcut, "time.scales" ) <- c( attr( Lcut, "time.scales" ), new.scales )
  attr( Lcut, "time.since"  ) <- c( attr( Lcut, "time.since"  ), new.states )
  }
    
# return the cut object without the auxilary variable
rmcol <- grep( "whseq", names(Lcut) )
Lcut[,-rmcol]  
  }






