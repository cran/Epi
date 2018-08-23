### R code from vignette source 'etm'
### Encoding: UTF-8

###################################################
### code chunk number 1: etm.rnw:42-52
###################################################
options( width=90 )
library( Epi )
library( etm )
library( survival )
library( splines )
print( sessionInfo(), l=F )
# load( url("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v038i04/dli.data.rda") )
# save( dli.data, file="./dli.Rda" )
load(           file="./dli.Rda" )
str( dli.data )


###################################################
### code chunk number 2: etm.rnw:55-56
###################################################
subset( dli.data, id %in% c(2,5,388,511,531,600) )


###################################################
### code chunk number 3: etm.rnw:61-62
###################################################
subset( dli.data, c(TRUE,diff(time)<0 & diff(id)==0 ) )


###################################################
### code chunk number 4: etm.rnw:65-66
###################################################
dli.data[dli.data$id==531 & dli.data$from==2,"time"] <- 0.45


###################################################
### code chunk number 5: etm.rnw:69-70
###################################################
with( dli.data, table( from, to ) )


###################################################
### code chunk number 6: etm.rnw:81-87
###################################################
dli <- transform( dli.data,
                  to = as.numeric(to),
                  ti = ave( time, id, FUN=function( x ) c(0,x[-length(x)]) ) )
dli$to <- ifelse( is.na(dli$to), dli$from, dli$to )
subset( dli, id %in% c(5,388,511,600) )
with( dli, table( from, to ) )


###################################################
### code chunk number 7: etm.rnw:102-108
###################################################
tmp <- subset(dli,to==2 & from!=to)[,c("id","time")]
names(tmp)[2] <- "tr"
dli <- merge( dli, tmp, all.x=TRUE )
dli$tr <- with( dli, ifelse( ti-tr>=0, ti-tr, NA ) )
subset(dli,id %in% c(600,603,608) )
str( dli )


###################################################
### code chunk number 8: etm.rnw:111-122
###################################################
# DLI
tmp <- subset(dli,to==4 & from!=to)[,c("id","time")]
names(tmp)[2] <- "tD"
dli <- merge( dli, tmp, all.x=TRUE )
dli$tD <- with( dli, ifelse( ti-tD>=0, ti-tD, NA ) )
# Rm2
tmp <- subset(dli,to==6 & from!=to)[,c("id","time")]
names(tmp)[2] <- "tR"
dli <- merge( dli, tmp, all.x=TRUE )
dli$tR <- with( dli, ifelse( ti-tR>=0, ti-tR, NA ) )
subset(dli,id %in% c(600,603,608) )


###################################################
### code chunk number 9: etm.rnw:129-142
###################################################
state.names <- c("Rem" , "D/Rem",
                 "Rel" , "D/Rel",
                 "DLI" , "D/DLI",
                 "Rem2", "D/Rem2",
                 "Rel2")
dli <- Lexis( entry        = list( tfi=ti, tfr=tr, tfD=tD, tfR=tR ),
              entry.status = factor( from, levels=0:8, labels=state.names ),
               exit        = list( tfi=time ),
               exit.status = factor(   to, levels=0:8, labels=state.names ),
                        id = id,
                      data = dli )
print.data.frame(
subset( dli, id %in% c(600,603,608) )[,1:13], digits=3)


###################################################
### code chunk number 10: etm.rnw:157-158
###################################################
summary( dli )


###################################################
### code chunk number 11: etm.rnw:164-165 (eval = FALSE)
###################################################
## boxes( dli )


###################################################
### code chunk number 12: boxes
###################################################
n.st <- nlevels( dli$lex.Cst )
# Colors for stages reflecting severity
st.col <- rep(c("limegreen","darkorange","yellow3","forestgreen","red"),each=2)[-10]
st.col[1:4*2] <- rgb( t(col2rgb(st.col[1:4*2])*0.5 + 255*0.5), max=255 )
boxes( dli, wmult=1.1, hmult=1.2, lwd=4,
            boxpos=list(x=c(10,30,30,50,50,70,70,90,90),
                        y=c(25, 8,42,25,59,42,76,59,93)),
            scale.R=100, show.BE=TRUE, DR.sep=c(" (",")"),
            col.bg=st.col, col.txt=rep(c("white","black"),5)[-10],
            col.border=c("white","black")[c(1,2,1,2,1,2,1,2,1)] )


###################################################
### code chunk number 13: stack
###################################################
st.dli <- stack( dli )
str( st.dli )
round(
cbind( "Original"=with(    dli, tapply( lex.dur, lex.Cst, sum ) ),
       "Stacked" =with( st.dli, tapply( lex.dur, lex.Cst, sum ) ) ), 1 )
round( xtabs( cbind( lex.Fail, lex.dur ) ~ lex.Tr, data = st.dli ), 1 )


###################################################
### code chunk number 14: Subset-mort
###################################################
dd.dli <- subset( st.dli, lex.Tr %in% levels(lex.Tr)[c(1,3,5,7)] )
table( dd.dli$lex.Tr )


###################################################
### code chunk number 15: etm.rnw:222-224
###################################################
dd.dli$lex.Tr <- factor( dd.dli$lex.Tr )
round( xtabs( cbind( lex.Fail, lex.dur ) ~ lex.Tr, data=dd.dli ), 1 )


###################################################
### code chunk number 16: Cox-0
###################################################
str( dd.dli )
c0 <- coxph( Surv(tfi,tfi+lex.dur,lex.Fail) ~ lex.Tr, data=dd.dli )
summary( c0 )


###################################################
### code chunk number 17: Cox-1
###################################################
dd.dli$Rst <- Relevel( dd.dli$lex.Tr, list(Remis=c(1,4),Relapse=2:3) )
with( dd.dli, table( lex.Tr, Rst ) )
c1 <- coxph( Surv(tfi,tfi+lex.dur,lex.Fail) ~ Rst, data=dd.dli )
summary( c1 )
anova( c0, c1, test="Chisq" )


###################################################
### code chunk number 18: split
###################################################
sp.dli <- splitLexis( dli, breaks=seq(0,50,1/10) )
print.data.frame( subset( sp.dli, id==603 )[,1:13], digits=3 )
summary( sp.dli )


###################################################
### code chunk number 19: stack-1
###################################################
ss.dli <- stack( sp.dli )
m.dli <- subset( ss.dli, lex.Tr %in% levels(lex.Tr)[c(1,3,5,7)] )
m.dli$lex.Tr <- factor( m.dli$lex.Tr )
m.dli$Rst <- Relevel( m.dli$lex.Tr, list(Remis=c(1,4),Relapse=2:3) )
with( m.dli, ftable( Rst, lex.Tr, lex.Fail, col.vars=3 ) )


###################################################
### code chunk number 20: tab-split-stack
###################################################
YDtab <- xtabs( cbind(lex.dur,lex.Fail) ~ I(floor(tfi*10)/10) + Rst,
                data=subset(m.dli,tfi<2.1) )
dnam <- dimnames(YDtab)
dnam[[3]] <- c("Y","D","rate")
YDrate <- array( NA, dimnames=dnam, dim=sapply(dnam,length) )
YDrate[,,1:2] <- YDtab
YDrate[,,3] <- YDrate[,,2]/YDrate[,,1]*1000
round( ftable( YDrate, row.vars=1 ), 1 )


###################################################
### code chunk number 21: m-knots
###################################################
( m.kn <- with( subset( m.dli, lex.Fail ),
                c(0,quantile( tfi+lex.dur, probs=1:4/5 )) ) )


###################################################
### code chunk number 22: Poisson-0
###################################################
m0 <- glm( lex.Fail ~ Ns( tfi, knots=m.kn ) + lex.Tr,
           family = poisson, offset=log(lex.dur), data=m.dli )
summary( m0 )
round(ci.lin( m0, E=T ),3)


###################################################
### code chunk number 23: Poisson-Cox-comp
###################################################
round(cbind( ci.lin( m0, subset="->" )[,1:2],
             ci.lin( c0              )[,1:2] ), 4 )
round(cbind( ci.lin( m0, subset="->" )[,1:2]/
             ci.lin( c0              )[,1:2] ), 4 )


###################################################
### code chunk number 24: Poisson-1
###################################################
m1 <- update( m0, . ~ . - lex.Tr + Rst )
anova( m0, m1, test="Chisq" )


###################################################
### code chunk number 25: time-interaction
###################################################
m0i <- update( m0, . ~ . + Ns( tfi, knots=m.kn ):lex.Tr )
m1i <- update( m1, . ~ . + Ns( tfi, knots=m.kn ):Rst )
anova( m0, m0i, test="Chisq" )
anova( m1, m1i, test="Chisq" )


###################################################
### code chunk number 26: Pairwise
###################################################
anova( m0i, m1i, test="Chisq" )


###################################################
### code chunk number 27: i-parms
###################################################
ci.lin( m1i )[,1:2]


###################################################
### code chunk number 28: pred-mort
###################################################
pr.pt <- seq(0,10,0.02)
 n.pt <- length( pr.pt )
CM <- cbind( 1, Ns( pr.pt, knots=m.kn ) )


###################################################
### code chunk number 29: mortality-rates
###################################################
Rem.Dead <- ci.exp( m1i, ctr.mat=cbind(CM,CM*0) )
Rel.Dead <- ci.exp( m1i, ctr.mat=cbind(CM,CM  ) )


###################################################
### code chunk number 30: rates
###################################################
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
matplot( pr.pt, cbind( Rem.Dead, Rel.Dead )*100, log="y", ylim=c(1,400),
         type="l", lty=1, lwd=c(3,1,1), las=1,
         col=rep(clr<-c("limegreen","red"),each=3),
         ylab="Mortality rates per 100 PY", xlab="Time since entry" )
text( c(10,10), 400*c(0.8,1), c("Remission","Relapse"), col=clr,
      adj=c(1,1) )


###################################################
### code chunk number 31: withRR
###################################################
RR.RmRl <- ci.exp( m1i, ctr.mat=cbind(CM*0,CM) )
par( mar=c(3,4,1,4), mgp=c(3,1,0)/1.6 )
matplot( pr.pt, cbind( Rem.Dead, Rel.Dead )*100, log="y", ylim=c(1,400),
         type="l", lty=1, lwd=c(3,1,1),
         col=rep(c("limegreen","red"),each=3),
         ylab="Rate per 100 PY", xlab="Time since entry", yaxt="n" )
abline( h=10 )
text( c(10,10), 400*c(0.8,1),
      c("Remission","Relapse"), col=c("limegreen","red"),
      adj=c(1,1), font=2 )
matlines( pr.pt, RR.RmRl*10, type="l", lty=1, lwd=c(3,1,1), col="blue" )
yt <- outer(c(1,2,5),c(1,10,100,1000),"*")
axis( side=2, at=yt, labels=yt   , las=1 )
axis( side=4, at=yt, labels=yt/10, las=1 )
mtext( "Rate ratio", side=4, line=2, col="blue" )


###################################################
### code chunk number 32: events
###################################################
str( ss.dli )
with( subset( ss.dli, lex.Fail & lex.Tr %in% levels(lex.Tr)[-c(1,3,5,7)] ),
      dotchart( tfi+lex.dur, groups=factor(lex.Tr), pch=16, cex=0.8) )


###################################################
### code chunk number 33: progress-subset
###################################################
p.dli <- subset( ss.dli, lex.Tr %in% levels(lex.Tr)[-c(1,3,5,7)] )
p.dli$lex.Tr <- factor( p.dli$lex.Tr )
( p.kn <-with( subset( p.dli, lex.Fail ),
               c(0,quantile( tfi+lex.dur, probs=1:3/4 )) ) )


###################################################
### code chunk number 34: etm.rnw:452-453
###################################################
p.dli$lex.Tr <- Relevel( p.dli$lex.Tr, list("Rem->Rel"=c(1,4), 2, 3 ) )


###################################################
### code chunk number 35: pr-modl
###################################################
p2 <- glm( lex.Fail ~ Ns( tfi, knots=p.kn, intercept=TRUE ):lex.Tr -1 ,
           family = poisson, offset=log(lex.dur), data=p.dli )
summary( p2 )
round(ci.exp( p2, Exp=FALSE ), 2 )


###################################################
### code chunk number 36: pr-ctr-mat
###################################################
CP <- Ns( pr.pt, knots=p.kn, intercept=TRUE )
Rem.Rel <- ci.exp( p2, subset="TrRem", ctr.mat=CP )
Rel.DLI <- ci.exp( p2, subset="TrRel", ctr.mat=CP )
DLI.Rem <- ci.exp( p2, subset="TrDLI", ctr.mat=CP )


###################################################
### code chunk number 37: pr-rates
###################################################
par( mar=c(3,4,1,4), mgp=c(3,1,0)/1.6 )
matplot( pr.pt, cbind(Rem.Rel,Rel.DLI,DLI.Rem)*100, log="y", ylim=c(1,500),
         type="l", lty=1, lwd=c(3,1,1),
         col=rep(c("limegreen","red","blue"),each=3),
         ylab="Progression rate per 100 PY", xlab="Time since entry" )
par(font=2)
legend( "topright", legend=c("Remission -> Relapse",
                             "Relapse -> DLI",
                             "DLI -> Remission"), bty="n",
        text.col=c("limegreen","red","blue"), xjust=1 )


###################################################
### code chunk number 38: pr-pt-interval
###################################################
pr.pt[1:5]
unique( diff(pr.pt) )
( il <- mean( diff(pr.pt) ) )


###################################################
### code chunk number 39: tr-mat
###################################################
states <- levels( dli$lex.Cst )
dnam <- list(From=states, To=states, time=pr.pt )
AR <- array( 0, dimnames=dnam, dim=sapply(dnam,length) )
str( AR )


###################################################
### code chunk number 40: fill-tr
###################################################
AR["Rem" ,"D/Rem" ,] <- Rem.Dead[,1]
AR["Rem2","D/Rem2",] <- Rem.Dead[,1]
AR["Rel" ,"D/Rel" ,] <- Rel.Dead[,1]
AR["DLI" ,"D/DLI" ,] <- Rel.Dead[,1]
AR["Rem" ,"Rel"   ,] <-  Rem.Rel[,1]
AR["Rem2","Rel2"  ,] <-  Rem.Rel[,1]
AR["Rel" ,"DLI"   ,] <-  Rel.DLI[,1]
AR["DLI" ,"Rem2"  ,] <-  DLI.Rem[,1]


###################################################
### code chunk number 41: diag-tr
###################################################
SI <- apply(AR,c(1,3),sum)


###################################################
### code chunk number 42: tp-mat
###################################################
AP <- AR
for( i in 1:(dim(AP)[3]) )
   {
   AP[,,i] <- AR[,,i]/SI[,i] * (1-exp(-SI[,i]*il))
   diag(AP[,,i]) <- exp(-SI[,i]*il)
   }
AP[is.na(AP)] <- 0
round( SI[,1], 4 )
round( ftable( AR[,,50+1:2], row.vars=c(3,1)), 4 )
round( ftable( AP[,,50+1:2], row.vars=c(3,1)), 4 )


###################################################
### code chunk number 43: state-occ
###################################################
pi0 <- c(1,rep(0,8))
ST <- SI*0
ST[,1] <- pi0 %*% AP[,,1]
for( i in 2:dim(ST)[2] ) ST[,i] <- ST[,i-1] %*% AP[,,i]
str( ST )
round(t(ST[,1:10]),3)


###################################################
### code chunk number 44: state-plot
###################################################
matplot( pr.pt, t(ST)[,-1], type="l", ylim=c(0,0.5), las=1,
         lty=1, lwd=2, col=rainbow(8),
         xlab="Time since 1st remission", ylab="Fraction of patients" )
legend( "topleft", legend=rownames(ST)[-1],
        bty="n", lty=1, lwd=2, col=rainbow(8),
        ncol=1 )


###################################################
### code chunk number 45: cum-states
###################################################
cST <- apply(ST,2,cumsum)


###################################################
### code chunk number 46: etm.rnw:612-613 (eval = FALSE)
###################################################
## matplot( pr.pt, t(cST), type="l", lty=1, lwd=2 )


###################################################
### code chunk number 47: state-occ-1
###################################################
perm <- c(1,3,5,7,9,2,4,6,8)
dimnames(ST)[[1]][perm]
matplot( pr.pt, t(apply(ST[perm,],2,cumsum)), type="l", lty=1, lwd=2 )


###################################################
### code chunk number 48: state-occ-2
###################################################
endpos <- cumsum(ST[perm,n.pt]) - ST[perm,n.pt]/2
matplot( pr.pt, t(apply(ST[perm,],2,cumsum)),
         type="l", lty=1, lwd=2, col=gray(c(0.6,0)[c(1,1,1,1,2,1,1,1,1)]),
         xlim=c(0,11.5), ylim=c(0,1), yaxs="i" )
text( 10.05, endpos, dimnames(ST)[[1]][perm], font=2, adj=0 )


###################################################
### code chunk number 49: state-occ-fun
###################################################
state.occ <-
function( perm=1:n.st, line=NULL )
{
clr <- st.col[perm]

mindist <- 1/40
endpos <- cumsum(ST[perm,n.pt]) - ST[perm,n.pt]/2
minpos <- (1:n.st-0.5)*mindist
maxpos <- 1-rev(minpos)
endpos <- pmin( pmax( endpos, minpos ), maxpos )

stkcrv <- t(apply(ST[perm,],2,cumsum))
matplot( pr.pt, stkcrv,
         type="l", lty=1, lwd=1, col="transparent",
         xlim=c(0,11.5), ylim=c(0,1), yaxs="i", xaxs="i", bty="n",
         xlab="Time since 1st remission", ylab="Fraction of patients" )
text( 10.05, endpos, dimnames(ST)[[1]][perm], font=2, adj=0, col=clr )
for( i in 9:1 )
polygon( c( pr.pt, rev(pr.pt) ),
         c( stkcrv[,i], if(i>1) rev(stkcrv[,i-1]) else rep(0,n.pt) ),
         col=clr[i], border=clr[i] )
if( !is.null(line) )
matlines( pr.pt, stkcrv[,line], type="l", lty=1, lwd=3, col="black" )
}


###################################################
### code chunk number 50: state-occ-col1
###################################################
state.occ( perm=1:n.st )


###################################################
### code chunk number 51: state-occ-col2
###################################################
state.occ( perm=c(1,3,5,7,9,2,4,6,8), line=5 )


###################################################
### code chunk number 52: etm.rnw:719-744
###################################################
get.rates <- function( N=10 )
{
Rem.Dead <- ci.lin( m1i,                 ctr.mat=cbind(CM,CM*0), sample=N )
Rel.Dead <- ci.lin( m1i,                 ctr.mat=cbind(CM,CM  ), sample=N )
Rem.Rel  <- ci.lin( p2 , subset="TrRem", ctr.mat=CP            , sample=N )
Rel.DLI  <- ci.lin( p2 , subset="TrRel", ctr.mat=CP            , sample=N )
DLI.Rem  <- ci.lin( p2 , subset="TrDLI", ctr.mat=CP            , sample=N )
states <- levels( dli$lex.Cst )
dnam <- list( From = states,
                To = states,
              time = pr.pt,
            sample = 1:N )
AR <- AP <- array( 0, dimnames=dnam, dim=sapply(dnam,length) )
AR["Rem" ,"D/Rem" ,,] <- exp(Rem.Dead)
AR["Rem2","D/Rem2",,] <- exp(Rem.Dead)
AR["Rel" ,"D/Rel" ,,] <- exp(Rel.Dead)
AR["DLI" ,"D/DLI" ,,] <- exp(Rel.Dead)
AR["Rem" ,"Rel"   ,,] <- exp(Rem.Rel )
AR["Rem2","Rel2"  ,,] <- exp(Rem.Rel )
AR["Rel" ,"DLI"   ,,] <- exp(Rel.DLI )
AR["DLI" ,"Rem2"  ,,] <- exp(DLI.Rem )
AR
}
system.time( AR <- get.rates(1000) )
str( AR )


###################################################
### code chunk number 53: AR-AP
###################################################
trans.prob <-
function( AR )
{
# A matrix for transition probabilities:
AP <- AR * 0
# Compute the interval length for the give rates
il <- mean( diff( as.numeric( dimnames(AR)[[3]] ) ) )
# Sum of the Intensities out of each state
SI <- apply(AR,c(1,3,4),sum)
for( i in 1:dim(AR)[3] ) # Loop over times
for( j in 1:dim(AR)[4] ) # Loop over samples
   {
   AP[,,i,j] <- AR[,,i,j]/SI[,i,j] * (1-exp(-SI[,i,j]*il))
   diag(AP[,,i,j]) <- exp(-SI[,i,j]*il)
   }
AP[is.na(AP)] <- 0
invisible( AP )
}
system.time( AP <- trans.prob( AR ) )


###################################################
### code chunk number 54: AP-ST
###################################################
pi0 <- rep(1:0,c(1,n.st-1))
ST  <- AP[1,,,]*0
names( dimnames(ST) )[1] <- "State"
str( ST )
system.time(
for( j in 1:dim(ST)[3] )
   {
   ST[,1,j] <- pi0 %*% AP[,,1,j]
   for( i in 2:dim(ST)[2] ) ST[,i,j] <- ST[,i-1,j] %*% AP[,,i,j]
   }
            )
str( ST )


###################################################
### code chunk number 55: def-st-oc-sim
###################################################
state.occ.sim <-
function( perm = 1:n.st,
           pct = c(5,95),
         cicol = rgb(1/9,1/9,1/9,1/9),
          line = NULL )
{
clr <- st.col[perm]
cST <- apply(ST[perm,,],2:3,cumsum)
cST <- apply(cST,1:2,quantile,probs=c(pct/100,0.5) )
mindist <- 1/40
endpos <- cST["50%",,n.pt] - diff(c(0,cST["50%",,n.pt]))/2
minpos <- (1:n.st-0.5)*mindist
maxpos <- 1-rev(minpos)
endpos <- pmin( pmax( endpos, minpos ), maxpos )

matplot( pr.pt, t(cST["50%",,]),
         type="n", # lty=1, lwd=2, col=gray(c(0.6,0)[c(1,2,1,2,1,2,1,2,1)]),
         xlim=c(0,11.5), ylim=c(0,1), yaxs="i", xaxs="i", bty="n",
         xlab="Time since 1st remission", ylab="Fraction of patients" )
text( 10.05, endpos, dimnames(ST)[[1]][perm], font=2, adj=0, col=clr )
for( i in n.st:1 )
polygon( c( pr.pt, rev(pr.pt) ),
         c( cST["50%",i,], if(i>1) rev(cST["50%",i-1,]) else rep(0,dim(cST)[3]) ),
         col=clr[i], border=clr[i] )
for( i in n.st:1 )
polygon( c( pr.pt, rev(pr.pt) ),
         c( cST[1,i,], rev(cST[2,i,]) ),
         col=cicol, border=cicol )
if( !is.null(line) ) matlines( pr.pt, t(cST["50%",c(line,NA),]),
                               type="l", lty=1, lwd=3, col="black" )
}


###################################################
### code chunk number 56: states-ci1
###################################################
par(mar=c(3,3,1,1),mgp=c(3,1,0)/1.6,las=1)
state.occ.sim( pct=c(5,95) )


###################################################
### code chunk number 57: states-ci2
###################################################
par(mar=c(3,3,1,1),mgp=c(3,1,0)/1.6,las=1)
state.occ.sim( perm=c(1,3,5,7,9,2,4,6,8), pct=c(5,95), line=5 )


###################################################
### code chunk number 58: states-cirem
###################################################
par(mar=c(3,3,1,1),mgp=c(3,1,0)/1.6,las=1)
state.occ.sim( perm=c(1,7,3,5,9,4,6,8,2), pct=c(5,95), line=c(2,5) )


###################################################
### code chunk number 59: parm-CLFS
###################################################
CLFS <- apply( ST["Rem",,]+ST["Rem2",,], 1, quantile, probs=c(500,25,975)/1000 )
str( CLFS )


###################################################
### code chunk number 60: CLFS
###################################################
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
matplot( pr.pt, t(CLFS),
         type="l", lty=1,lwd=c(3,1,1), col="black",
         ylim=c(0,1),
         ylab="P(CLFS)", xlab="Time since first remission" )


###################################################
### code chunk number 61: etm-ex
###################################################
tra <- matrix(FALSE, 9, 9,
              dimnames = list(as.character(0:8), as.character(0:8)))
tra[1, 2:3] <- TRUE
tra[3, 4:5] <- TRUE
tra[5, 6:7] <- TRUE
tra[7, 8:9] <- TRUE
### computation of the transition probabilities
dli.etm <- etm::etm(dli.data, as.character(0:8), tra, "cens", s = 0)
str(dli.etm)
### Computation of the clfs + var clfs
clfs <- dli.etm$est["0", "0", ] +
        dli.etm$est["0", "6", ]
var.clfs <- dli.etm$cov["0 0", "0 0", ] +
            dli.etm$cov["0 6", "0 6", ] +
        2 * dli.etm$cov["0 0", "0 6", ]
## computation of the 95% CIs + plot
ciplus  <- clfs + qnorm(0.975) * sqrt(var.clfs)
cimoins <- clfs - qnorm(0.975) * sqrt(var.clfs)

plot(dli.etm$time, clfs, type = "s", lwd=3,
     bty = "n", ylim = c(0, 1), yaxs="i", las=1,
     xlab = "Time since 1st remission (years)",
     ylab = "P(CLFS)" )
lines(dli.etm$time, cimoins, lty = 3, type = "s")
lines(dli.etm$time, ciplus , lty = 3, type = "s")
matlines( pr.pt, t(CLFS), lty=c(1,3,3), lwd=c(3,1,1), col="red" )


###################################################
### code chunk number 62: etm-ex2
###################################################
xdli.etm <- etm( dli )
str( xdli.etm )
plot( xdli.etm, col=rainbow(15), lty=1, lwd=3,
      legend.pos="topright", bty="n", yaxs="i" )


