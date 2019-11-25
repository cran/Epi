list.lex <-
function( Lx, dig, sh="_" )
    {
            wh <- c("lex.id",timeScales(Lx),"lex.dur","lex.Cst","lex.Xst")
           oth <- setdiff(names(Lx),wh)
         tmpLx <- Lx[,c(wh,oth)]
       wh.num  <- sapply(tmpLx,is.numeric)
tmpLx[,wh.num] <- round(tmpLx[,wh.num],dig)
  names(tmpLx) <- gsub("lex.",sh,names(tmpLx))
print( tmpLx, row.names=FALSE )
return( NULL )
    }
