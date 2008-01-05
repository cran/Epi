doCutLexis <- function(data, cut, timescale) {
    
    ## Code each new interval using new variable lex.cut:
    ## 0 = unchanged interval (cut occurs after exit)
    ## 1 = first part of split interval
    ## 2 = second part of split interval (or cut occurs before entry)
    
    cut[is.na(cut)] <- Inf #If a cut time is missing, it never happens

    ## First intervals (before the cut)
    in.1 <- entry(data, timescale)
    ex.1 <- pmin(cut, exit(data, timescale))

    ## Create Lexis object for first intervals
    lx.1 <- data
    lx.1$lex.dur  <- ex.1 - in.1
    lx.1$lex.cut <- ifelse(cut < exit(data, timescale), 1, 0)
    
    ## Second intervals (after the cut)
    in.2 <- pmax(cut, entry(data, timescale))
    ex.2 <- exit(data, timescale)

    ## Create Lexis object for second intervals
    lx.2 <- data
    lx.2$lex.dur <- ex.2 - in.2
    lx.2$lex.cut <- 2
    
    ## Update entry times
    lx.2[, timeScales(data)] <- exit(data) - lx.2$lex.dur
    
    return(rbind(lx.1, lx.2))
}

setStatus.default <- function(data, new.state)
{
    data$lex.Xst[data$lex.cut == 1] <- new.state[data$lex.cut == 1]
    data$lex.Cst[data$lex.cut == 2] <- new.state
    
    return(data)
}

setStatus.numeric <- function(data, new.state, precursor.states=NULL,
                              progressive=TRUE) {
    
    if (!is.numeric(new.state)) {
        stop("If lex.Cst, lex.Xst are numeric, new.state must be numeric too")
    }
    
    data$lex.Xst[data$lex.cut == 1] <- new.state[data$lex.cut == 1]
    data$lex.Cst[data$lex.cut == 2] <- new.state

    exit.state <- data$lex.Xst[data$lex.cut == 2]
    is.precursor <- exit.state %in% precursor.states
    if (progressive) {
        is.precursor <- is.precursor | (exit.state < new.state)
    }
    data$lex.Xst[data$lex.cut == 2][is.precursor] <- new.state[is.precursor]

    return(data)
}

setStatus.factor <- function(data, new.state, precursor.states=NULL,
                             progressive=TRUE) {
    
    if (!is.character(new.state)) {
        stop("new.state must be a character vector, but it is ",str(new.state))
    }
    
    current.states <- levels(data$lex.Cst)
    new.states <- setdiff(new.state, current.states)

    ## Modify factor levels if necessary
    if (length(new.states) > 0) {
        all.states <- c(current.states, sort(new.states))
        levels(data$lex.Cst) <- all.states
        levels(data$lex.Xst) <- all.states
    }

    data$lex.Xst[data$lex.cut == 1] <- new.state[data$lex.cut == 1]
    data$lex.Cst[data$lex.cut == 2] <- new.state

    exit.state <- data$lex.Xst[data$lex.cut==2]
    is.precursor <-  exit.state %in% precursor.states
    if (progressive) {
        if (is.ordered(data$lex.Xst)) {
            is.precursor <- is.precursor | (exit.state < new.state)
        }
        else {
            warning("progressive=TRUE argument ignored for unordered factor")
        }
    }
    data$lex.Xst[data$lex.cut==2][is.precursor] <- new.state[is.precursor]

    return(data)
}

# Added by BxC
match.cut <-
function( data, cut )
   {
   if( sum(!is.na(match(c("lex.id","cut","new.state"),names(cut))))<3 )
       stop( "The dataframe supplied in the cut= argument must have columns",
             "'lex.id','cut','new.state', but the columns are:\n", names( cut ) )
   else
     {
     if( length( unique( cut$lex.id ) ) < nrow( cut ) )
         stop( "Values of 'lex.id' must be unique in the 'cut' dataframe" )
     else
       zz <- merge( data[,"lex.id",drop=FALSE], cut, all.x=TRUE )
       if( is.factor ( data$lex.Cst ) ) zz$new.state <- as.character(zz$new.state)
       if( is.numeric( data$lex.Cst ) ) zz$new.state <- as.numeric(zz$new.state)
       return( zz )
     }
   }
# End of addition / change


cutLexis <- function(data,
                     cut,
                     timescale = 1,
                     new.state,
                     progressive = FALSE,
                     precursor.states = NULL,
                     count = FALSE)

{
    if (!inherits(data, "Lexis"))
      stop("data must be a Lexis object")

# Added by BxC: Enabling the count argument in cutLexis.
    if( count )
      return( countLexis( data=data, cut=cut, timescale=timescale ) )
# End of addition

# Added by BxC
    if( inherits( cut, "data.frame" ) ){
      zz <- match.cut( data, cut )
      cut <- zz$cut
      new.state <- zz$new.state
    }
    else if (length(cut) == 1) {
# End of addition / change
     cut <- rep(cut, nrow(data))
    }
    else if (length(cut) != nrow(data)) {
        stop("'cut' must have length 1 or nrow(data) (=", nrow(data),
             "),\n --- but it has length ", length(cut),".")
    }

    timescale <- Epi:::check.time.scale(data, timescale)
    if (length(timescale) > 1) {
        stop("Multiple time scales")
    }

    if (missing(new.state)) {
        new.state <- data$lex.Cst       #Carry forward last state
    }
    else if (length(new.state) == 1) {
        new.state <- rep(new.state, nrow(data))
    }
    else if (length(new.state) != nrow(data)) {
        stop("'new.state' must have length 1 or nrow(data) (=", nrow(data),
             "),\n --- but it has length ", length(new.state))
    }

    if (progressive) {
        if (is.factor(data$lex.Cst) && !is.ordered(data$lex.Cst)) {
            stop("progressive=TRUE invalid for unordered factors")
        }
        if (any(data$lex.Xst < data$lex.Cst)) {
            stop("Lexis object is not progressive before splitting")
        }
    }
    
    lx <- doCutLexis(data, cut, timescale)
    if (is.factor(data$lex.Cst)) {
        lx <- setStatus.factor(lx, new.state, precursor.states, progressive)
    }
    else if (is.numeric(data$lex.Cst)) {
        lx <- setStatus.numeric(lx, new.state, precursor.states, progressive)
    }
    else {
        lx <- setStatus.default(lx, new.state)
    }
    
    ## Remove redundant intervals
    lx <- lx[lx$lex.dur > 0,]

# BxC addition
    ## Remove the lex.cut column
    lx <- lx[,-match("lex.cut",names(lx))]
# end addition

    attributes( lx ) <- attributes( data )
    return(lx[order(lx$lex.id,lx[,timescale]),])
}


countLexis <- function(data, cut, timescale = 1)
{
    if (!inherits(data, "Lexis"))
      stop("data must be a Lexis object")

# Added by BxC
    if( inherits( cut, "data.frame" ) ){
      zz <- match.cut( data, cut )
      cut <- zz$cut
      new.state <- zz$new.state
    }
    else if (length(cut) == 1) {
# End of addition / change
        cut <- rep(cut, nrow(data))
    }
    else if (length(cut) != nrow(data)) {
        stop("'cut' must have length 1 or nrow(data) (=", nrow(data),
             "),\n --- but it has length ", length(cut),".")
    }

    timescale <- Epi:::check.time.scale(data, timescale)
    if (length(timescale) > 1) {
        stop("Multiple time scales")
    }
    
    lx <- doCutLexis(data, cut, timescale)

    ## Update status variables
    lx$lex.Xst[lx$lex.cut == 1] <- lx$lex.Cst[lx$lex.cut == 1] + 1
    lx$lex.Cst[lx$lex.cut == 2] <- lx$lex.Cst[lx$lex.cut == 2] + 1
    lx$lex.Xst[lx$lex.cut == 2] <- lx$lex.Xst[lx$lex.cut == 2] + 1
    
    ## Remove redundant intervals
    lx <- lx[lx$lex.dur > 0,]

# BxC addition
    ## Remove the lex.cut column
    lx <- lx[,-match("lex.cut",names(lx))]
# end addition

    return(lx[order(lx$lex.id,lx[,timescale]),])
}
