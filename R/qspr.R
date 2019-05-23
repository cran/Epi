spread <-
function( dd, tt, rand=TRUE )
{
dd <- tapply(dd,tt,sum,na.rm=TRUE)
tt <- as.numeric(names(dd))
oo <- order(tt)
dd <- dd[oo]
tt <- tt[oo]
# now d is the number of items scored with value t
dt <- diff(tt)
nt <- length(tt)
# here are the nt+1 boundaries of the nt intervala we shall use
tx <- tt[c(1,1:nt)] + c(-dt[1],dt[c(1:(nt-1),(nt-1))])/2
# the vector to hold the spread out values of dd
dx <- numeric(sum(dd))
cd <- c(0,cumsum(dd))
for( i in 1:nt ) dx[(cd[i]+1):cd[i+1]] <-
                 if( rand ) runif(dd[i],tx[i],tx[i+1])
                 else seq(tx[i],tx[i+1],,nn<-dd[i]+2)[-c(1,nn)]
return( sort(dx) )
}
dd <- c(5,4,6,2,10)
tt <- c(1,3,4,7,11)

( dr <- spread(dd,tt) )
( ds <- spread(dd,tt,rand=FALSE) )
plot(tt,dd,type="h",lwd=4,xlim=c(0,15),ylim=c(0,max(dd)))
abline(v=quantile(rep(tt,dd),1:5/6),lty=3)
points( ds, rep(1.5,length(ds)), pch=16, col="red" )
abline(v=quantile(ds,1:5/6),col="red")
points( dr, rep(0.5,length(ds)), pch=16, col="blue" )
abline(v=quantile(dr,1:5/6),col="blue")
