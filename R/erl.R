erl <-
function( int,
          muW,
          muD,
          lam = NULL,
       age.in = 0,
            A = NULL,
       immune = is.null(lam),
          yll = TRUE,
         note = TRUE )
{
# Computes expected residual life time for Well and Dis states
# respectively in an illness-death model, optionally ignoring
# the well->ill transition

# Utility to integrate a survival function from the last point where
# it is 1, assuming points are 1 apart
trsum <-
function( x )
{
x[c(diff(x)==0,TRUE)] <- NA     
sum( ( x[-length(x)] + x[-1] ) / 2, na.rm=TRUE )
}

# Check sensibility
if( !immune & is.null(lam) ) stop( "'lam' is required when immune=FALSE\n" )

# Survival functions    
             sD <- surv1( int=int,      muD,      age.in = age.in, A = A )
if( immune ) sW <- surv1( int=int, muW,           age.in = age.in, A = A )    
else         sW <- surv2( int=int, muW, muD, lam, age.in = age.in, A = A )

# Area under the survival functions    
erl <- cbind( apply( sW[,-1], 2, trsum ),
              apply( sD[,-1], 2, trsum ) ) * int
colnames( erl ) <- c("Well","Dis")
rownames( erl ) <- colnames( sW )[-1]

# Should we compute years of life lost? 
if( yll ) erl <- cbind( erl, YLL = erl[,"Well"] - erl[,"Dis"] )

# Cautionary note
if( immune )
  {
  attr( erl, "NOTE" ) <- "Calculations assume that Well persons cannot get Ill (quite silly!)."
  if( note ) cat("NOTE:",  attr( erl, "NOTE" ), "\n" )
  }  
return( erl )
}

# yll is just a simple wrapper for erl, only selecting the YLL column.
yll <-
function( int,
          muW,
          muD,
          lam = NULL,
       age.in = 0,
            A = NULL,
       immune = is.null(lam),
         note = TRUE ) erl( int = int,      
                            muW = muW,      
                            muD = muD,      
                            lam = lam,      
                         age.in = age.in,   
                              A = A,        
                         immune = immune,   
                            yll = TRUE,
                           note = note )[,"YLL"]

surv1 <-
function( int, mu, age.in=0, A=NULL )
{
# Computes the survival function from age A till the end, assuming
# that mu is a vector of mortalities in intervals of length int. 
# int and mu should be in compatible units that is T and T^-1 for
# some unit T (months, years, ...) 

# age-class boundaries    
age <- 0:length(mu)*int + age.in

# cumulative rates and survival at the boundaries    
Mu <- c( 0, cumsum( mu )*int )
Sv <- exp( -Mu )
surv <- data.frame( age=age, surv=Sv )
    
# if a vector of conditioning ages A is given
if( cond <- !is.null(A) )
  {
  j <- 0
  # actual conditioning ages
  cage <- NULL
  for( ia in A )
     {
     j <- j+1
     # Where is the age we condition on
     cA <- which( diff(age>ia)==1 )
     surv <- cbind( surv, pmin( 1, surv$surv/(surv$surv[cA]) ) )
     cage[j] <- surv$age[cA]
     }
  }
names( surv )[-1] <- paste( "A", c( age.in, if( cond ) cage else NULL ), sep="" )
rownames( surv ) <- NULL
return( surv ) 
}

erl1 <-
function( int, mu, age.in = 0 )
{
# Computes expected residual life time at all ages    
age <- 0:length(mu)*int + age.in

# Small utility: cumulative cumulative sum from the end of a vector
musmuc <- function( x ) rev( cumsum( rev(x) ) )

# The survival function with a 0 at end, and the integral from the upper end
surv <- surv1( int = int, mu = mu, age.in = age.in )[,2]
cbind( age = age,
      surv = surv,
       erl = c( musmuc( ( surv[-1]-diff(surv)/2 ) ) /
                surv[-length(surv)], 0 ) * int )
}

surv2 <-
function( int, muW, muD, lam, age.in=0, A=NULL )
{
# check the vectors
if( length(muW) != length(muD) | 
    length(muD) != length(lam) ) 
  stop( "Vectors with rates must have same length:\n",
          "length(muW)=", length(muW),
        ", length(muD)=", length(muD),
        ", length(lam)=", length(lam) )

# First the workhorse that computes the survival function for a
# person in Well assuming that the mortality rate from this state is
# muW, disease incidence is in lam, and mortality in the diseased
# state is muD, and that all refer to constant rates intervals of
# length int starting from age.in, conditional on survival to A 
wsurv2 <-
function( int, muW, muD, lam, age.in=0, A=0 )
{
# age-class boundaries - note one longer that rate vectors refers to
# boundaries of intervals not midpoints  
age <- 0:length(muW)*int + age.in
    
# cumulative rates at the boundaries, given survival to A
MuW <-  cumsum( c( 0, muW ) * ( age > A ) ) * int
MuD <-  cumsum( c( 0, muD ) * ( age > A ) ) * int
Lam <-  cumsum( c( 0, lam ) * ( age > A ) ) * int
    
# probability of being well
pW <- exp( -( Lam + MuW ) )
    
# probability of diagnosis at s --- first term in the integral for
# P(DM at a). Note that we explicitly add a 0 at the start so we get a
# probability of 0 of transition at the first age point 
Dis <- c(0,lam) * ( age > A ) * exp( -(Lam+MuW) ) * int
    
# for each age (age[ia]) we compute the integral over the range
# [0,age] of the product of the probability of diagnosis and the
# probability of surviving from diagnosis till age ia
pDM <- Dis * 0
for( ia in 1:length(age) )
   pDM[ia] <- sum( Dis[1:ia] * exp( -(MuD[ia]-MuD[1:ia]) ) )
                   # 1st term as function of s (1:ia)
                               # 2nd term integral over range s:age
                   # upper integration limit is age (ia) and the lower
                   # limit is the intermediate age (at DM) (1:ia)
# Finally, we add the probabilities of being in Well resp. DM to get
# the overall survival:
surv <- data.frame( age = age, surv = pDM + pW )
return( surv )
}    

# survival from start    
surv <- wsurv2( int, muW, muD, lam, age.in=age.in, A=0 )

# add columns for conditioning ages    
if( !is.null(A) )
   { 
   for( j in 1:length(A) )
      { 
      surv <- cbind( surv, 
                    wsurv2( int, muW, muD, lam, age.in=age.in, A=A[j] )[,2] )
      }
   }
Al <- A
for( i in 1:length(A) ) Al[i] <- max( surv$age[surv$age <= A[i]] )
colnames( surv )[-1] <- paste( "A", c( age.in, Al ), sep="" )    

# done!    
return( surv )
}
