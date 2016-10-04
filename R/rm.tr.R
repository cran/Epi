rm.tr <-
function( obj, from, to )
    {    
    # checks
    if( from==to) stop( "Don't be silly, 'from' and 'to' are identical." )
    if( !(from %in% levels(obj) ) ) stop( "'from' not a state in the object." )
    if( !(  to %in% levels(obj) ) ) stop( "  'to' not a state in the object." )
### These things do not change over the purging iterations 
    # Sort the rows of the object and count them   
    obj <- obj[order(obj$lex.id,obj[,timeScales(obj)[1]]),]
    nrw <- nrow( obj )
    # First obs for each person    
    no1 <- !duplicated( obj$lex.id )
    wh1 <- which( no1 )
### Utility function doing the work inside the loop        
purge.one <-
function( obj, from, to, nrw, wh1 )
    {
    # Where are the illegal transitions 
    chX <- ( paste( obj$lex.Cst, obj$lex.Xst ) == paste( from, to ) )
    whX <- which( chX )
    if( length(whX)>0 )
      {  
      # Change lex.Xst in this record
      obj$lex.Xst[whX] <- obj$lex.Cst[whX]
      # and lex.Cst in next record, but only if it is not a new person   
      whX <- setdiff( whX[whX<nrw]+1, wh1 ) 
      obj$lex.Cst[whX] <- obj$lex.Cst[whX-1]
      return( obj )
      }
    else return(NULL)
    }
### Iterate till all illegal transitions are weeded out        
cont <- TRUE
while( cont )
    {
    zzz <- obj
    cont <- !is.null( obj <- purge.one( obj=zzz,
                                        from=from, to=to,
                                        nrw = nrw, wh1 = wh1 ) )
    }
### Return last actual object        
zzz
}

