tab <-
function (x, ...)
UseMethod("tab")

tab.Lexis <-
function( x, simplify=TRUE, scale=1, ... )
{
tr <- trans <- with( x, table(lex.Cst,lex.Xst) )
for( i in intersect(rownames(trans),colnames(trans)) ) tr[i,i] <- 0
trans <- addmargins(trans)
tr    <- addmargins(tr)
tr    <- tr[,ncol(tr)]
pyrs  <- with( x, addmargins( tapply(lex.dur,lex.Cst,sum,na.rm=TRUE),
                              FUN=function(x) sum(x,na.rm=TRUE) ) )
res <- cbind( trans, tr, pyrs, tr/pyrs, tr/pyrs/exp(1.96/sqrt(tr)),
                                        tr/pyrs*exp(1.96/sqrt(tr)) )
res[,ncol(res)-2:0] <- res[,ncol(res)-2:0]*scale
colnames( res )[ncol(res)-4:0] <- c(" #events:"," #risk time:",
                             "Rate", "   (95%", " c.i.)")
names( dimnames( res ) ) <- c("From","\nStates:\n     #records:\n     To")
if( simplify ) res <- res[!is.na(pyrs),]
if( nrow( res )==2 ) res <- res[1,,drop=FALSE]
res
}
