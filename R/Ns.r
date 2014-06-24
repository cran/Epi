#-------------------------------------------------------------------------------
# A wrapper for ns that automatically  takes the smallest and largest knots
# as boundary knots without further ado
Ns <-
function( x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots=NULL )
{
 if( !is.null(df) ) ns( x, df=df )
 else
   {
   if( is.null( Boundary.knots ) )
     {
     if( !is.null( knots ) )
       {
       knots <- sort( unique( knots ) )
       ok <- c(1,length(knots))
       Boundary.knots <- knots[ok]
       knots <- knots[-ok]
       }
     }
   ns( x, knots=knots, intercept=intercept, Boundary.knots=Boundary.knots )
   }
}

