Lexis.lines <-
function( entry.date = NA, 
           exit.date = NA,
          birth.date = NA,
           entry.age = NA,
            exit.age = NA,
           risk.time = NA, 
            col.life = "black",
            lwd.life = 2,
                fail = NA,
            cex.fail = 1.1,
            pch.fail = c(1, 16),
            col.fail = rep(col.life, 2),
                data = NULL
          ) 
{
if( !is.null( data ) )
  {
  attach( data, 2 )
  on.exit( detach( pos=2 ) )
  }

# Complete the information on lifelines
XX <- Life.lines( entry.date = entry.date,
                   entry.age = entry.age,
                   exit.date = exit.date,
                    exit.age = exit.age,
                   risk.time = risk.time,
                  birth.date = birth.date )
  
# Was XX returned as a Date-object?
# If so make a numerical version i LL, otherwise just a copy.
#
if( attr( XX, "Date" ) )
  {
  LL <- data.frame( lapply( XX, unclass ) )
  LL[,c(1,3,5)] <- LL[,c(1,3,5)] / 365.25 + 1970
  LL[,c(2,4,6)] <- LL[,c(2,4,6)] / 365.25
  } else LL <- XX

# Find age and date ranges in the current plot.
#
date <- par( "usr" )[1:2]
age  <- par( "usr" )[3:4]

# Truncate life-lines to fit inside the plot.
# delay is the risk time at start falling outside (dealyed entry)
# early is the risk time at end falling outside (early retirement)
#
delay <- pmax( 0, date[1]-LL[,1], age[1]-LL[,2] )
early <- pmax( 0, LL[,3]-date[2], LL[,4]-age[2] )

# Remove those falling entirely outside the diagram
#
out <- ( delay + early ) > LL[,6]

# If there are any left at all we plot them inside the diagram
#
if( sum( out ) < length( out ) )
  {
  segments( (LL[,1] + delay)[!out], 
            (LL[,2] + delay)[!out], 
            (LL[,3] - early)[!out], 
            (LL[,4] - early)[!out], lwd=lwd.life, col=col.life )
  if( is.numeric( fail ) ) fail <- ( fail > 0 )
  points( LL[,3], LL[,4],
          pch=16,
          col=par()$bg,
          cex=cex.fail )
  points( LL[,3], LL[,4],
          pch=pch.fail[(fail & early==0)+1],
          col=col.fail[(fail & early==0)+1],
          cex=cex.fail )
  }

# Return the untouched version of the completed dataframe
#
invisible( data.frame( XX, fail=fail ) )
}
