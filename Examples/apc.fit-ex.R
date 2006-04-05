library( Epi )
data(lungDK)

# Taylor a dataframe that meets the requirements
exd <- lungDK[,c("Ax","Px","D","Y")]
names(exd)[1:2] <- c("A","P")

# Two different ways of parametrizing the APC-model, ML
ex.H <- apc.fit( exd, npar=7, model="ns", drift="Holford",  parm="ACP", scale=10^5 )
ex.W <- apc.fit( exd, npar=7, model="ns", drift="weighted", parm="ACP", scale=10^5 )

# Sequential fit, first AC, then P given AC.
ex.S <- apc.fit( exd, npar=7, model="ns", parm="AC-P", scale=10^5 )

# Show the estimated drifts
ex.H[["Drift"]]
ex.W[["Drift"]]
ex.S[["Drift"]]

# First nice plot frame
par( mar=c(3,4,1,4), mgp=c(3,1,0)/1.5, las=1 )
sc <- apc.frame(a.lab = seq( 30, 90, 20 ),
                a.tic = seq( 30, 90, 10 ),
               cp.lab = seq( 1860, 2000, 20 ),
               cp.tic = seq( 1860, 2000, 10 ),
                r.lab = c(1,2,5,10,20,50)*10,
                r.tic = c(1:9*10, 1:5*100),
               rr.ref = 200,
                  gap = 22 )

# Reference lines
abline( v=ex.H[[5]][2]-sc[1] )
segments( 1860-sc[1], sc[2], 2000-sc[1], sc[2] )

# Fill in the estimated effects
apc.lines( ex.S, col="green", scale="rates", frame.par=sc, ci=TRUE )
apc.lines( ex.H, col="red", scale="rates", frame.par=sc, ci=TRUE )
apc.lines( ex.W, col="blue", scale="rates", frame.par=sc, ci=TRUE )

# Extract the drifts in \% per year
rd <- formatC( (ex.W[["Drift"]]["A-d",]-1)*100, format="f", digits=2 )
wd <- formatC( (ex.W[["Drift"]]["APC",]-1)*100, format="f", digits=2 )
hd <- formatC( (ex.H[["Drift"]]["APC",]-1)*100, format="f", digits=2 )

# Put them on the plot
text( 1999-sc[1], 11,
      paste( "Weighted drift:", wd[1], "(", wd[2], "-", wd[3], ") \%/year" ),
      col="blue", adj=c(1,0), font=2, cex=0.8 )
text( 1999-sc[1], 11*1.2,
      paste( "Naïve drift:", hd[1], "(", hd[2], "-", hd[3], ") \%/year" ),
      col="red", adj=c(1,0), font=2, cex=0.8 )
text( 1999-sc[1], 11*1.2^2,
      paste( "Raw drift:", rd[1], "(", rd[2], "-", rd[3], ") \%/year" ),
      col="green", adj=c(1,0), font=2, cex=0.8 )

