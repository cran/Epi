library( Epi )
data( nickel )

stat.table(  cut( agein, breaks=c(0,50,60,75) ),
   contents=list( N=count(),
                  Y=sum( ageout-agein ),
                  D=sum( icd %in% c(162,163) ) ) )

stat.table(  cut( agein, breaks=c(0,50,60,75) ),
   contents=list( N=count(),
                  Y=sum( ageout-agein ),
                 mb=mean( dob ),
                  D=sum( icd %in% c(162,163) ) ) )

stat.table(  cut( agein, breaks=c(0,50,60,75) ),
   contents=list( N=count(),
                  Y=sum( ageout-agein ),
                 mb=mean( dob ),
               Rate=ratio( icd %in% c(162,163), ageout-agein, 1000 ) ) )

stat.table(  cut( agein, breaks=c(0,50,60,75) ),
   contents=list( N=count(),
                  Y=sum( ageout-agein ),
               Rate=ratio( icd %in% c(162,163), ageout-agein, 1000 ) ) )

