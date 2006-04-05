# Workhorse 3:
# ------------
# Function to cut the follow-up from "enter" to "exit" at point "dof"
# and put a failure indicator at that point.
fcut1 <-
function( enter,
           exit,
           fail,
            dof,
     fail.value = 1,
           data = data.frame( enter, exit, fail, dof ),
         Expand = 1:nrow( data ) )
{
dorg <- data
attach( dorg )
on.exit( detach( dorg ) )
dof[is.na(dof)] <- Inf

# Dataframe of first intervals (new failures)
Enter <- enter
Exit  <- pmin( dof, exit, na.rm=T )
d1 <- data.frame( Expand=Expand, Enter=Enter, Exit=Exit,
                                 Fail=ifelse( Exit==dof, fail.value, fail ), data )
d1 <- d1[Enter<Exit & enter<dof,]

# Dataframe of last intervals (possibly failures already known)
Enter <- pmax( dof, enter, na.rm=TRUE )
Exit  <- exit
d2 <- data.frame( Expand=Expand, Enter=Enter, Exit=Exit,
                                 Fail=fail.value, data )
d2 <- d2[Enter<Exit & exit>dof,]

# Combine the two and sort them
dx <- rbind( d1, d2 )
dx[order(dx$Expand),]
}

