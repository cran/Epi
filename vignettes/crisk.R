### R code from vignette source 'crisk'
### Encoding: UTF-8

###################################################
### code chunk number 1: crisk.rnw:24-27
###################################################
options( width=90,
         SweaveHooks=list( fig=function()
         par(mar=c(3,3,1,1),mgp=c(3,1,0)/1.6,las=1,bty="n") ) )


###################################################
### code chunk number 2: crisk.rnw:83-100
###################################################
library(Epi)
library(popEpi)
data(DMlate)
Ldm <- Lexis(entry = list( per = dodm,
                           age = dodm-dobth, 
                           tfd = 0 ),
              exit = list( per = dox ),
       exit.status = factor( !is.na(dodth), labels = c("DM","Dead") ),
              data = DMlate )
summary(Ldm, t = T)
Mdm <- mcutLexis( Ldm,
                   wh = c('dooad','doins'),
           new.states = c('OAD','Ins'),
            precursor = 'DM',
           seq.states = FALSE,
                 ties = TRUE )
summary( Mdm )


###################################################
### code chunk number 3: crisk.rnw:105-107
###################################################
Sdm <- splitMulti(factorize(subset(Mdm, lex.Cst == "DM")),
                  tfd = seq(0, 20, 1/12))


###################################################
### code chunk number 4: boxes5
###################################################
boxes(Mdm, boxpos = list(x = c(15, 50, 15, 85, 85),
                         y = c(85, 50, 15, 85, 15)), 
          scale.R = 100, 
          show.BE = TRUE)


###################################################
### code chunk number 5: boxes4
###################################################
boxes( Relevel(Sdm, c(1, 4, 2, 3)), 
       boxpos  = list(x = c(15, 85, 80, 15),
                      y = c(85, 85, 20, 15)),
       scale.R = 100, 
       show.BE = TRUE )


###################################################
### code chunk number 6: crisk.rnw:139-142
###################################################
mD <- gam.Lexis(Sdm, ~ s(tfd, k = 5), to = 'Dead')
mO <- gam.Lexis(Sdm, ~ s(tfd, k = 5), to = 'OAD' )
mI <- gam.Lexis(Sdm, ~ s(tfd, k = 5), to = 'Ins' )


###################################################
### code chunk number 7: crisk.rnw:155-160
###################################################
int <- 1/100
nd <- data.frame( tfd = seq(int,10,int)-int/2 ) # not the same as the split, 
                                                # and totally unrelated to it
rownames(nd) <- nd$tfd
str(nd)


###################################################
### code chunk number 8: rates
###################################################
matshade(nd$tfd, cbind(ci.pred(mD, nd),
                       ci.pred(mI, nd),
                       ci.pred(mO, nd))*1000, 
         ylim = c(0.02,500), yaxt = "n",
         ylab = "Rates per 1000 PY", 
         xlab = "Time since DM diagnosis (years)",
         col = c("black","red","blue"), log = "y", lwd = 3, plot = TRUE)
axis(side = 2, at = ll<-outer(c(1,2,5),-2:3,function(x,y) x*10^y),
               labels = formatC(ll,digits = 4), las = 1)
axis(side = 2, at = ll<-outer(c(1.5,2:9),-2:3,function(x,y) x*10^y),
               labels = NA, tcl = -0.3)
text(0, 0.5*0.6^c(1,2,0), 
     c("Dead","Ins","OAD"),
     col = c("black","red","blue"), adj = 0)


###################################################
### code chunk number 9: crisk.rnw:194-209
###################################################
# rates at midpoints
lD <- ci.pred( mD, nd )[,1]
lI <- ci.pred( mI, nd )[,1]
lO <- ci.pred( mO, nd )[,1]
# cumulative rates and survival fuction at right border of the intervals
LD <- cumsum(lD) * int
LI <- cumsum(lI) * int
LO <- cumsum(lO) * int
Sv <- exp( -LD - LI - LO )
# but when integrating to get the cumulative risks we use the average
# of the survival function at the two endpoints (adding 1 as the first)
mp <- function(x) x - diff(c(1, x)) / 2
rD <- cumsum(lD * mp(Sv)) * int
rI <- cumsum(lI * mp(Sv)) * int
rO <- cumsum(lO * mp(Sv)) * int


###################################################
### code chunk number 10: crisk.rnw:214-218
###################################################
summary(rD + rI +rO + Sv)
oo <- options(digits = 20)
cbind(summary(Sv + rD + rI + rO))
options(oo)


###################################################
### code chunk number 11: stack
###################################################
zz <- mat2pol(cbind(rD,rI,rO,Sv), x = nd$tfd, 
              xlim = c(0,10), xaxs = "i", yaxs = "i", las = 1,
              xlab = "Time since DM diagnosis (years)", 
              ylab = "Probability",
               col =  c("black","red","blue","forestgreen") )
mm <- t(apply(zz,1,mid<-function(x) x[-1]-diff(x)/2))
text( 9, mm[900,], c("Dead","Ins","OAD","DM"), col = "white" )
box(col = "white",lwd = 3)


###################################################
### code chunk number 12: crisk.rnw:278-280
###################################################
head(cbind(ci.pred(mI,nd),     ci.exp(mI,nd)            ))
head(cbind(ci.pred(mI,nd), exp(ci.lin(mI,nd)[,c(1,5:6)])))


###################################################
### code chunk number 13: crisk.rnw:285-287
###################################################
str(ci.lin(mI, nd, sample = 4))
head(cbind(ci.pred(mI,nd), exp(ci.lin(mI, nd, sample = 4))))


###################################################
### code chunk number 14: crisk.rnw:319-321 (eval = FALSE)
###################################################
## setwd("/home/bendix/stat/R/examples")
## source('ci.Crisk.R', echo=TRUE, max=10000)


###################################################
### code chunk number 15: crisk.rnw:326-334
###################################################
system.time(
res <- ci.Crisk(list(OAD = mO, 
                     Ins = mI, 
                    Dead = mD),
                            nd = data.frame(tfd = (1:1000-0.5)/100),
                            nB = 1000,
                          perm = 4:1))
str(res)


###################################################
### code chunk number 16: crisk.rnw:364-372
###################################################
system.time(
rsm <- ci.Crisk(list(OAD = mO, 
                     Ins = mI, 
                    Dead = mD),
                            nd = data.frame(tfd = (1:1000-0.5)/100),
                            nB = 2000,
                       sim.res = 'rates'))
str(rsm) 


###################################################
### code chunk number 17: crisk.rnw:379-387
###################################################
system.time(
csm <- ci.Crisk(list(OAD = mO, 
                     Ins = mI, 
                    Dead = mD),
                            nd = data.frame(tfd = (1:1000-0.5)/100),
                            nB = 2000,
                       sim.res = 'crisk'))
str(csm) 


###################################################
### code chunk number 18: crisk.rnw:406-408
###################################################
Brates <- aperm(apply(rsm, 1:2, Epi:::mnqt), c(2,3,1))
str(Brates)


###################################################
### code chunk number 19: rates-ci
###################################################
matshade(nd$tfd, cbind(ci.pred(mD, nd),
                       ci.pred(mI, nd),
                       ci.pred(mO, nd))*1000, 
         ylim = c(0.1,500), yaxt = "n",
         ylab = "Rates per 1000 PY", 
         xlab = "Time since DM diagnosis (years)",
         col = c("black","red","blue"), log = "y", lwd = 3, plot = TRUE)
matlines(nd$tfd, cbind(Brates[,"Dead",],
                       Brates[,"Ins" ,],
                       Brates[,"OAD" ,])*1000,
         col = c("white","black","black"), lty = 3, lwd=c(3,1,1))
axis(side = 2, at = ll<-outer(c(1,2,5),-2:3,function(x,y) x*10^y),
               labels = formatC(ll,digits = 4), las = 1)
axis(side = 2, at = ll<-outer(c(1.5,2:9),-2:3,function(x,y) x*10^y),
               labels = NA, tcl = -0.3)
text(0, 0.5*0.6^c(1,2,0), 
     c("Dead","Ins","OAD"),
     col = c("black","red","blue"), adj = 0)


###################################################
### code chunk number 20: crates
###################################################
matshade(c(0,nd$tfd+1/200),
         cbind(res$Crisk[,"Dead",],
               res$Crisk[,"Ins" ,],
               res$Crisk[,"OAD" ,]), plot = TRUE,
         xlim = c(0,10), xaxs = "i", yaxs = "i", las = 1,
         xlab = "Time since DM diagnosis (years)", 
         ylab = "Cumulative probability",
          col = c("black","red","blue"))
text(8, 0.3 + c(1,0,2)/25, 
     c("Dead","Ins","OAD"),
     col = c("black","red","blue"), adj = 0)


###################################################
### code chunk number 21: crisk.rnw:467-469
###################################################
str(res$Crisk)
str(res$Srisk)


###################################################
### code chunk number 22: stack-ci
###################################################
zz <- mat2pol(res$Crisk[,c("Dead","Ins","OAD","Surv"),1],
              x = as.numeric(dimnames(res$Crisk)[[1]])/100,
           xlim = c(0,10), xaxs = "i", yaxs = "i", las = 1,
           xlab = "Time since DM diagnosis (years)", 
           ylab = "Probability",
            col =  c("black","red","blue","forestgreen") )
mm <- t(apply(zz, 1, mid<-function(x) x[-1] - diff(x) / 2))
text( 9, mm[900,], c("Dead","Ins","OAD","DM"), col = "white" )
matshade(as.numeric(dimnames(res$Srisk)[[1]])/100,
         cbind(res$Srisk[,1,],
               res$Srisk[,2,],
               res$Srisk[,3,]),
         col = 'transparent', col.shade = "white", alpha = 0.3)


###################################################
### code chunk number 23: crisk.rnw:503-504
###################################################
str(res$Stime)


###################################################
### code chunk number 24: crisk.rnw:507-510
###################################################
s510 <- res$Stime[1:2*500,,]
dimnames(s510)[[1]] <- c(" 5 yr","10 yr")
round(ftable(s510, row.vars=1:2), 2)


