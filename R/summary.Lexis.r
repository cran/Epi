summary.Lexis <-
function( object, simplify=TRUE, scale=1, by=NULL,
          Rates=FALSE, timeScales=FALSE, ... )
{
# If we have a by argument find out what to do
    if (!is.null(by)) {
        if (is.character(by)) {
            if (!all(by %in% names(object)))
                stop("Wrong 'by' argument: '", paste(by,collapse="','"),
                     "' - must be name(s) of variable(s) in the Lexis object")
            else res <- lapply(split(object, object[, by]), summary.Lexis,
                               by = NULL, simplify = simplify, scale = scale,
                               Rates = Rates, timeScales=timeScales, ...)
        }
        else {
            if (length(by) != nrow(object))
                stop("Wrong length of 'by' argument:", length(by),
                     "must be same length as rows of the Lexis object:", nrow(object) )
            else res <- lapply(split(object, by), summary.Lexis,
                               by = NULL, simplify = simplify, scale = scale,
                               Rates = Rates, timeScales=timeScales, ...)
        }
        # to avoid printing the time scale information repeatedly
        for( i in 1:(length(res)-1) ) res[[i]]$timeScales <- NULL
        return( res )
    }

# Table(s) of all transitions (no. records)
tr <- trans <- with( object, table(lex.Cst,lex.Xst) )

# Remove diagonal, i.e. records with no transition
for( i in intersect(rownames(trans),colnames(trans)) ) tr[i,i] <- 0

# Margins added
trans <- addmargins(trans)
tr    <- addmargins(tr)    # Sum omitting the diagonal
trm   <- tr[,ncol(tr)]

# Compute person-years in each Cst-state
pyrs  <- with( object,
               addmargins( tapply( lex.dur,lex.Cst,sum,na.rm=TRUE),
                                   FUN=function(x) sum(x,na.rm=TRUE) ) )/scale
pers  <- with( object,
               c( tapply( lex.id, lex.Cst, function(x) length(unique(x)) ),
                  length(unique(lex.id)) ) )

# Amend the table of records with column of events and person-years
trans <- cbind( trans, trm, pyrs, pers )

# Annotate the table nicely
colnames( trans )[ncol(trans)-2:0] <-
    c(" Events:","Risk time:"," Persons:" )
colnames( trans )[ncol(tr)] <- " Records:"
names( dimnames( trans ) ) <- c("From","\nTransitions:\n     To")

# Make the rates and annotate the table nicely
rates <- sweep( tr, 1, pyrs, "/" )
colnames( rates )[ncol(rates)] <- "Total"
names( dimnames( rates ) ) <-
     c("From",
       paste("\nRates",
             if( scale != 1 ) paste(" (per ",scale,")",sep=""),
             ":\n     To", sep="") )
if( simplify )
  {
  trans <- trans[!is.na(pyrs),]
  rates <- rates[!is.na(pyrs),]
  }
if( nrow(trans)==2 )
  trans <- trans[1,,drop = FALSE]
res <- list( Transitions = trans,
                   Rates = rates[-nrow(rates),,drop=FALSE],
              timeScales = data.frame( "time scale" = attr( object, "time.scales" ),
                                       "time since" = attr( object,  "time.since" ) )
)
if( !timeScales ) res <- res[-3]
if( !Rates      ) res <- res[-2]
class( res ) <- "summary.Lexis"
res
}

print.summary.Lexis <-
function( x, ..., digits=2 )
{
print( round( x$Transitions, digits ) )
if( "Rates" %in% names(x) )
print( round( x$Rates      , digits ) )
# if( "timeScales" %in% names(x) )
if( !is.null(x$timeScales) )
  {  
cat("\nTimescales:\n")
print( x$timeScales )
  }
}
