rcutLexis <-
function( Lx,
         cut,
   timescale = 1,
precursor.states )
{
# avoid note about no visible binding
new.state <- NULL
   lex.id <- NULL

# All new states shoudl be precursor states
pr.st <- unique(c(precursor.states, unique(cut$new.state)))
# utility to slect n'th element in vector
getn <- function(x,n) x[n]
# make the data frame a sorted, grouped tibble
# note that the cut in arrage is the variable cut in the tibble
cut <- group_by(cut, lex.id)
cut <- arrange(cut, cut, .by_group = TRUE)
# the resulting Lexis object 
nxL <- Lx
# max no transitions for one person
maxn <- max(table(cut$lex.id))
# loop over transition number per person
for(n in 1:maxn)
   {
   # the n'th transitions for each person
   nx.cut <- summarize( cut,
                        cut = getn(      cut, n),
                  new.state = getn(new.state, n) )
   # update with cuts at these for each person
   nxL <- cutLexis( nxL,  
                    cut = nx.cut,
              timescale = timescale,
       precursor.states = pr.st )
   }
# done, return result
nxL
}
