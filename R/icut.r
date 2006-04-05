# Workhorse 2:
# ------------
# Function to cut the follow-up from "en" to "ex" at point "cut"
# and carry fail on to the last.
icut <-
function( enter,
           exit,
            cut,
           fail = 0,
     cens.value = 0,
           data = data.frame( enter, exit, fail, cut ),
         Expand = 1:nrow( data ),
         na.cut = Inf )
{
# Attach a copy of the dataframe
dorg <- data
attach( dorg )
on.exit( detach( dorg ) )

# Missing cutpoints. Default is to treat then as +Inf
cut[is.na(cut)] <- na.cut

# Dataframe of first intervals
Enter <- enter
Exit  <- pmin( cut, exit, na.rm=T )
d1 <- data.frame( Expand=Expand, Enter=Enter, Exit=Exit,
                                 Fail=ifelse(exit<=Exit, fail, cens.value ),
                                 Time=0, dorg )
d1 <- d1[Enter<Exit & Enter<cut,]
# Dataframe of last intervals
Enter <- pmax( cut, enter, na.rm=TRUE )
Exit  <- exit
d2 <- data.frame( Expand=Expand, Enter=Enter, Exit=Exit,
                                 Fail=fail, Time=1, dorg )
d2 <- d2[Enter<Exit & exit>cut,]

dx <- rbind( d1, d2 )
dx[order(dx$Expand,dx$Time),]
}

