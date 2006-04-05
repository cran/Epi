# Workhorse 5:
# ------------
# Function to expand follow-up along one time-scale
ex1 <-
function( enter,
           exit,
           fail,
         origin = 0,
          scale = 1,
         breaks,
           data = data.frame( enter, exit, fail ),
         Expand = 1:nrow(data) )
{
# Function to expand along one timescale

# Attach a copy of the dataframe if given
dorg <- data.frame( Expand=Expand, data )
attach( dorg )
on.exit( detach( dorg ) )

# 0-row dataframe for the results
dres <- data.frame( Expand=NA, Enter=NA, Exit=NA, Fail=NA, Time=NA )[NULL,]

# Loop through the intervals using isec
for( i in 2:length(breaks) )
   {
   ints <- isec( enter = (enter-origin)/scale,
                  exit = (exit-origin)/scale,
                  fail = fail,
                   int = breaks[(i-1):i],
                Expand = Expand )
      if( nrow( ints ) > 0 )
   dres <- rbind( dres, data.frame( ints, Time=breaks[i-1] ) )
   }
res <- merge( dres, dorg, by="Expand" )
res[order(res$Expand,res$Time),]
}

