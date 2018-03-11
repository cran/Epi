######################################################################
bsLexis <-
function( Tr, # List of lists of transition objects
        init, # Lexis object of persons to simulate.
           N = 1, # No. persons simulated per line in init
          Nb = 10, # No. of bootstrap samples
      lex.id,
     t.range = 20, # Range for rate computation in the simulation
       n.int = 101, # length of time intervals
    time.pts = seq(0,t.range,length.out=n.int)
         )
{
lapply( Tr, names )    
}
