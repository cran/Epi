Pplot <-
function( rates,
          age = as.numeric( dimnames( rates )[[1]] ),
          per = as.numeric( dimnames( rates )[[2]] ),
         grid = FALSE,
       p.grid = grid,
        ygrid = grid,
     col.grid = gray( 0.9 ),
        p.lim = range( per, na.rm=TRUE ) + c(0,diff(range(per))/30),
         ylim = range( rates, na.rm=TRUE ),
        p.lab = names( dimnames( rates ) )[2],
         ylab = deparse( substitute( rates ) ),
           at = NULL,
       labels = paste( at ),
         type = "l",
          lwd = 2,
          lty = 1,
          col = par( "fg" ),
       log.ax = "y",
          las = 1,
          ann = FALSE,
      cex.ann = 0.8,
        xannx = 1/20,
       a.thin = seq( 1, length( age ), 2 ),
          ...
          )
{
# Plot the frame
if( ann ) p.lim <- p.lim + c(0,diff( range( age ) ) * xannx )
matplot( per, t(rates), type="n",
         xlim=p.lim, ylim=ylim, xlab=p.lab, ylab=ylab,
         log=log.ax, las=las, if( !is.null( at ) ) yaxt="n" )
if( !is.null( at ) ) axis( side=2, at=at, labels=labels )
# and the grid if required
if( !missing( p.grid ) | !missing( grid ) )
  {
  if( is.logical( p.grid ) ) p.grid <- nice( per, log=par("xlog") )
  abline( v=p.grid, col=col.grid )
  }
if( !missing( ygrid ) | !missing( ygrid ) )
  {
  if( is.logical( ygrid ) ) ygrid <- nice( rates[!is.na(rates)], log=par("ylog") )
  abline( h=ygrid, col=col.grid )
  }
box()
# then the curves
matlines( per, t(rates), lwd=lwd, lty=lty, col=col, type=type, ... )
# annotate them if required (every second by default )
if( ann )
  {
   nr <- nrow( rates )
   nc <- ncol( rates )
   text( rep( per[nc], nr )[a.thin], rates[,nc][a.thin],
         paste( "", age[a.thin] ), adj=c(0,0.5),
         cex=cex.ann, col=if( length(col)==1 ) col else col[a.thin] ) 
  }
}
