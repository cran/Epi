pctab <- function( TT, margin=length( dim( TT ) ) )
  {
  nd <- length( dim( TT ) )
  sw <- (1:nd)[-margin[1]]
  sweep( addmargins( TT,
                     margin,
                     list( list( All=sum,
                                   N=function( x ) sum( x )^2/100 ) ) ),
         sw,
         apply( TT,
                sw,
                sum )/100,
         "/" )
  }
