summary.Lexis <-
function( object, simplify=TRUE, scale=1, ... )
{
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
res <- list( Transitions=trans, Rates=rates[-nrow(rates),,drop=FALSE] )
class( res ) <- "summary.Lexis"
res
}

print.summary.Lexis <-
function( x, ..., digits=2 )
{
print( round( x$Transitions, digits ) )
print( round( x$Rates      , digits ) )
}
