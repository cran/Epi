# Workhorse 4:
# ------------
# Function to cut the follow-up from "en" to "ex" at multiple
# points given in the list "dof".
# Basically loops over the event times and calls fcut1 each time.
fcut <-
function( enter,
           exit,
            dof,
           fail = 0,
           data = data.frame( enter, exit ),
         Expand = 1:nrow( data ) )
{
dorg <- data
attach( dorg )
on.exit( detach( dorg ) )

if( !(is.numeric( dof ) | is.list( dof ) ) )
  stop( "dof must be either a numerical vector or a list of vectors" )
  
if( is.numeric( dof ) ) dof <- list( dof )

# Dataframe with the columns of failure times at the end, and
# a variable counting the number of events so far at the start
dx <- data.frame( Expand=Expand,
                  Enter=enter, Exit=exit, Fail=fail, n.Fail=0, dof )
nf <- length( dof )
# Loop over the list of failure time vectors.
for( j in 1:nf )
   {
   # The variables Expand, Enter, Exit and Fail are formed anew at
   # each loop iteration by fcut1, and new values input for the next.
   # The purpose of using the data-argument is to have
   # fcut1 expand the later event times
   omit <- match( c("Expand","Enter","Exit","Fail"),
                  names( dx ) )
   dx <- fcut1( enter = dx$Enter,
                 exit = dx$Exit,
                 fail = dx$Fail,
                  dof = dx[,ncol(dx)-nf+j],
               Expand = dx$Expand,
                 data = dx[,-omit] )
   }
# A counter of no. events so far; if we have an interval with
# time after event and Enter equal to the currently
# processed event time
dx$n.Fail <- apply( dx[,(ncol(dx)-nf+1):ncol(dx),drop=F] <= dx$Enter,
                    1, sum, na.rm=TRUE )
# Remove the columns with failure times added
res <- dx[,1:(ncol(dx)-nf)]
# Merge the dataframe onto the expanded dataset
dx <- data.frame( Expand=Expand, dorg )
merge( res, dx, by="Expand", sort=TRUE )
}

