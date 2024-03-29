rcutLexis <-
function(Lx, cut,
       timescale = 1,
precursor.states = transient(Lx))
{
# avoid note about no visible binding
new.state <- lex.id <- NULL

# All new states should be precursor states
# new.state may be factor, hence the as.character
pr.st <- unique(c(precursor.states,
                  unique(as.character(cut$new.state))))

# utility to select n'th element in a vector
getn <- function(x,n) x[n]
    
# make the data frame a sorted, grouped tibble
# note that the cut in arrange is the variable cut in the tibble
cut <- group_by(cut, lex.id) %>%
       arrange(cut, .by_group = TRUE)
    
# nxL is the Lexis object to return eventually 
nxL <- Lx
    
# max no transitions for any one person in cut
maxn <- max(table(cut$lex.id))
    
# loop over transition number per person
for(n in 1:maxn)
   {
   # the n'th transitions for each person
   nx.cut <- summarize(cut,
                       cut = getn(      cut, n),
                 new.state = getn(new.state, n))
   # update with cuts at these for each person
   nxL <- cutLexis(nxL,  
                   cut = nx.cut,
             timescale = timescale,
      precursor.states = pr.st )
   }
# done, return result
nxL
}
