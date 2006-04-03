# Workhorse 1:
# ------------
# Function to determine the intersection between
# the intervals (enter,exit) and the FIXED interval int
isec <-
function( enter,
           exit,
           fail = 0,
            int,
     cens.value = 0,
         Expand = 1:length( enter ) )
{
 Enter <- pmax( int[1], enter )
 Exit  <- pmin( int[2], exit )
 ## The ">=" instead of "==" excludes overlap when enter>exit
 overlap <- ( Enter < Exit )
 Fail <- ifelse( exit > int[1] & exit <= int[2], fail, cens.value )
 # Collect result in matrix
 res <- data.frame( Expand=Expand, Enter=Enter, Exit=Exit, Fail=Fail )
 res[overlap,,drop=FALSE]
}

