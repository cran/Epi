cutLexis <-
function( data,
           cut,
     timescale = timeScales(data)[1],
     new.state = max( as.integer( c(factor(data$lex.Cst),factor(data$lex.Xst)) )+1 ),
        na.cut = Inf,
          cens = ifelse( is.character(new.state), "0", 0 ),
         count = missing(new.state) & missing(cens) )
{
force( count )
if( !inherits( data, "Lexis" ) )
    stop("First argument must be a Lexis object, but it is a ", class(data) )
if( !( timescale %in% (tscal<-timeScales(data)) ) )
    stop("You must specify one of the valid timescales: ", tscal )
    
# If a dataframe is supplied, match to get the entire vector of cuts
if( inherits( cut, "data.frame" ) )
  {
  if( "lex.id" %in% names( cut ) )
    {
    if( any( wh <- match( names(cut), timeScales(data) ) ) )
      {
      # Tease out the timescale, and change name to "cut"
      wh <- wh[!is.na(wh)]
      if( length(wh)>1 ) stop("More than one matching timescale in 'cut'")
      timescale <- timeScales( data )[wh]
      names( cut ) <- gsub( timescale, "cut", names(cut) )
      }
    if( "cut" %in% names( cut ) )
      {
      wh <- match( data$lex.id, cut$lex.id )
      new.state <- cut$new.state[wh]
      cut <- cut$cut[wh]
      }
    else
    stop( "'cut' or one of the timescales must be among the columns in the\n",
          "dataframe supplied in the cut= argument, they are:\n",
          names( cut ) )
    }
  else
  stop( "'lex.id' must be among the columns in the\n",
        "dataframe supplied in the cut= argument. The columns are:\n",
        names( cut ) )
  }
else
  {
  if( length( cut )==1 ) cut <- rep(cut,nrow(data))
  if( length( cut ) != nrow( data ) )
      stop( "'cut' must have length 1 or nrow(data) (=", nrow(data),
            "),\n --- but it has length ", length(cut),"." )
  }
if( length( new.state )==1 ) new.state <- rep(new.state,nrow(data))
if( length( new.state ) != nrow( data ) )
    stop( "'new.state' must have length 1 or nrow(data) (=", nrow(data),
          "),\n --- but it has length ", length(new.state) )
if( count & (!is.numeric(data$lex.Cst)|!is.numeric(data$lex.Xst)) )
    stop( "State variables lex.Cst and lex.Xst must be numeric" )

if( !count )
{
# Transform all to factors with the same set of levels
if( !is.character( new.state ) ) new.state <- as.character( new.state )
new.states <- unique( new.state )
data$lex.Cst <- factor( data$lex.Cst )
data$lex.Xst <- factor( data$lex.Xst )
all.states <- unique(c(levels(data$lex.Cst),
                       levels(data$lex.Xst),
                       new.states))

# Then convert the factors to numeric for the calculations
     new.state <- match(    new.state, all.states )
data$lex.Cst   <- match( data$lex.Cst, all.states )
data$lex.Xst   <- match( data$lex.Xst, all.states )
cens <- match( as.character(cens), all.states )
}
else cens <- sort(unique( c(data$lex.Cst, data$lex.Xst) ) )

# Missing cutpoints. Default is to treat as +Inf
cut[is.na(cut)] <- na.cut

# First intervals
in.1 <- entry( data, timescale )
ex.1 <- pmin( cut, as.vector( exit( data, timescale ) ), na.rm=T )
lx.1 <- data
lx.1$lex.dur  <- as.vector(ex.1 - in.1)
if( !count )
{    
lx.1$lex.Xst <- ifelse( ex.1 == as.vector( exit( data, timescale ) ),
                        data$lex.Xst,
                             new.state )
}
else
{    
lx.1$lex.Xst <- ifelse( ex.1 == as.vector( exit( data, timescale ) ),
                        data$lex.Xst,
                        data$lex.Cst+1 )
}

# Last intervals
in.2 <- pmax( cut, as.vector( entry( data, timescale ) ), na.rm=TRUE )
ex.2 <- exit( data, timescale )
lx.2 <- data
# Length of 2nd intervals i.e. AFTER the cut
lx.2$lex.dur  <- as.vector(ex.2 - in.2)
# Update timescale starting points by subtracting from the endpoint
lx.2[,timeScales(lx.2)] <- exit(data) - lx.2$lex.dur

# Note: we don't have to bother about intervals entirely to the left of cut,
# they will have dur < 0 and be removed below.

if( !count )
{    
# New entry status: if the interval contains a cut or
#                   if entire interval starts in a censoring state
lx.2$lex.Cst <- ifelse( ( lx.2$lex.dur < data$lex.dur | lx.2$lex.Cst %in% cens ),
                             new.state,
                        data$lex.Cst )
# Exit status is the new state if the original record had no transition from
# a censoring state
lx.2$lex.Xst <- ifelse( data$lex.Xst %in% cens,
                             new.state,
                        data$lex.Xst )
}
else
{    
# New entry status: if the interval contains a cut or
#                   if entire interval starts after the cut
lx.2$lex.Cst <- ifelse( ( lx.2$lex.dur < data$lex.dur | lx.2$lex.Cst %in% cens ),
                        data$lex.Cst+1,
                        data$lex.Cst )
# Exit status is the new state if the original record had no transition from
# a censoring state
lx.2$lex.Xst <- ifelse( data$lex.Xst %in% cens,
                        data$lex.Cst+2,
                        data$lex.Xst )
}
    
# Combine and restrict to the records with positive follow-up
lx <- rbind( lx.1, lx.2 )
lx <- lx[lx$lex.dur>0,]
if( !count )
{    
lx$lex.Cst <- factor( lx$lex.Cst, levels=1:length(all.states), labels=all.states )
lx$lex.Xst <- factor( lx$lex.Xst, levels=1:length(all.states), labels=all.states )
}
return( lx[order(lx$lex.id,lx[,timescale]),] )
}
