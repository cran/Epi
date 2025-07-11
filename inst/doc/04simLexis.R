### R code from vignette source '04simLexis.rnw'

###################################################
### code chunk number 1: 04simLexis.rnw:29-41
###################################################
options(width = 90,
        show.signif.stars = FALSE,
        SweaveHooks=list(fig = function()
                         par(mar = c(3, 3, 1, 1),
                             mgp = c(3, 1, 0) / 1.6,
                             las = 1,
                            lend = "butt",
                             bty = "n")))
library(Epi)
library(popEpi)
library(survival)
clear()


###################################################
### code chunk number 2: 04simLexis.rnw:44-46
###################################################
anfang <- Sys.time()
cat("Start time:", format(anfang, "%F, %T"), "\n")


###################################################
### code chunk number 3: 04simLexis.rnw:48-54
###################################################
vers <-
data.frame(R = substr(R.version.string, 11, 15),
         Epi = as.character(packageVersion(   "Epi")),
      popEpi = as.character(packageVersion("popEpi")))
names(vers) <- paste(" ", names(vers))
print(vers, row.names = FALSE)


###################################################
### code chunk number 4: start
###################################################
options( width=90 )
library( Epi )
print( sessionInfo(), l=F )


###################################################
### code chunk number 5: Lexis
###################################################
data(DMlate)
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit = list(Per=dox),
        exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
               data = DMlate )


###################################################
### code chunk number 6: cut
###################################################
dmi <- cutLexis( dml, cut = dml$doins,
                      pre = "DM",
                new.state = "Ins",
                new.scale = "t.Ins",
             split.states = TRUE )
summary( dmi, timeScales=T )


###################################################
### code chunk number 7: boxes
###################################################
boxes( dmi, boxpos = list(x=c(20,20,80,80),
                        y=c(80,20,80,20)),
            scale.R = 1000, show.BE = TRUE )


###################################################
### code chunk number 8: split
###################################################
Si <- splitLexis( dmi, seq(0,20,1/4), "DMdur" )
summary( Si )
print( subset( Si, lex.id==97 )[,1:10], digits=6 )


###################################################
### code chunk number 9: knots
###################################################
nk <- 5
( ai.kn <- with( subset(Si,lex.Xst=="Ins" & lex.Cst!=lex.Xst ),
                 quantile( Age+lex.dur  , probs=(1:nk-0.5)/nk ) ) )
( ad.kn <- with( subset(Si,lex.Xst=="Dead"),
                 quantile( Age+lex.dur  , probs=(1:nk-0.5)/nk ) ) )
( di.kn <- with( subset(Si,lex.Xst=="Ins" & lex.Cst!=lex.Xst ),
                 c(0,quantile( DMdur+lex.dur, probs=(1:(nk-1))/nk ) )) )
( dd.kn <- with( subset(Si,lex.Xst=="Dead"),
                 c(0,quantile( DMdur+lex.dur, probs=(1:(nk-1))/nk ) )) )
( ti.kn <- with( subset(Si,lex.Xst=="Dead(Ins)"),
                 c(0,quantile( t.Ins+lex.dur, probs=(1:(nk-1))/nk ) )) )


###################################################
### code chunk number 10: Poisson
###################################################
library( splines )
DM.Ins <- glm( (lex.Xst=="Ins") ~ Ns( Age  , knots=ai.kn ) +
                                  Ns( DMdur, knots=di.kn ) +
                                  I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
ci.exp( DM.Ins )
class( DM.Ins )


###################################################
### code chunk number 11: 04simLexis.rnw:308-314
###################################################
DM.Ins <- glm.Lexis( Si, from = "DM", to = "Ins",
                      formula = ~ Ns( Age  , knots=ai.kn ) +
                                  Ns( DMdur, knots=di.kn ) +
                                  I(Per-2000) + sex )
ci.exp( DM.Ins )
class( DM.Ins )


###################################################
### code chunk number 12: 04simLexis.rnw:319-328
###################################################
DM.Dead <- glm.Lexis( Si, from = "DM", to = "Dead",
                       formula = ~ Ns( Age  , knots=ad.kn ) +
                                   Ns( DMdur, knots=dd.kn ) +
                                   I(Per-2000) + sex )
Ins.Dead <- glm.Lexis( Si, from = "Ins",
                        formula = ~ Ns( Age  , knots=ad.kn ) +
                                    Ns( DMdur, knots=dd.kn ) +
                                    Ns( t.Ins, knots=ti.kn ) +
                                    I(Per-2000) + sex )


###################################################
### code chunk number 13: prop-haz
###################################################
All.Dead <- glm.Lexis( Si, to = c("Dead(Ins)","Dead"),
                      formula = ~ Ns( Age  , knots=ad.kn ) +
                                  Ns( DMdur, knots=dd.kn ) +
                                  lex.Cst +
                                  I(Per-2000) + sex )
round( ci.exp( All.Dead ), 3 )


###################################################
### code chunk number 14: get-dev
###################################################
what <- c("null.deviance","df.null","deviance","df.residual")
( rD <- unlist(  DM.Dead[what] ) )
( rI <- unlist( Ins.Dead[what] ) )
( rA <- unlist( All.Dead[what] ) )
round( c( dd <- rA-(rI+rD), "pVal"=1-pchisq(dd[3],dd[4]+1) ), 3 )


###################################################
### code chunk number 15: pr-array
###################################################
pr.rates <- NArray( list( DMdur = seq(0,12,0.1),
                          DMage = 4:7*10,
                          r.Ins = c(NA,0,2,5),
                          model = c("DM/Ins","All"),
                           what = c("rate","lo","hi") ) )
str( pr.rates )


###################################################
### code chunk number 16: mknd
###################################################
nd <- data.frame( DMdur = as.numeric( dimnames(pr.rates)[[1]] ),
                lex.Cst = factor( 1, levels=1:4,
                                  labels=levels(Si$lex.Cst) ),
                    sex = factor( 1, levels=1:2, labels=c("M","F")) )


###################################################
### code chunk number 17: make-pred
###################################################
for( ia in dimnames(pr.rates)[[2]] )
   {
dnew <- transform( nd, Age = as.numeric(ia)+DMdur,
                       Per = 1998+DMdur )
pr.rates[,ia,1,"DM/Ins",] <- ci.pred(  DM.Dead, newdata = dnew )
pr.rates[,ia,1,"All"   ,] <- ci.pred( All.Dead, newdata = dnew )
for( ii in dimnames(pr.rates)[[3]][-1] )
   {
dnew = transform( dnew, lex.Cst = factor( 2, levels=1:4,
                                          labels=levels(Si$lex.Cst) ),
                          t.Ins = ifelse( (DMdur-as.numeric(ii)) >= 0,
                                           DMdur-as.numeric(ii), NA ) )
pr.rates[,ia, ii ,"DM/Ins",] <- ci.pred( Ins.Dead, newdata = dnew )
pr.rates[,ia, ii ,"All"   ,] <- ci.pred( All.Dead, newdata = dnew )
    }
    }


###################################################
### code chunk number 18: mort-int
###################################################
par( mar=c(3,3,1,1), mgp=c(3,1,0)/1.6, las=1 )
plot( NA, xlim=c(40,82), ylim=c(5,300), bty="n",
      log="y", xlab="Age", ylab="Mortality rate per 1000 PY" )
abline( v=seq(40,80,5), h=outer(1:9,10^(0:2),"*"), col=gray(0.8) )
for( aa in 4:7*10 ) for( ii in 1:4 )
   matshade( aa+as.numeric(dimnames(pr.rates)[[1]]),
             cbind( pr.rates[,paste(aa),ii,"DM/Ins",],
                    pr.rates[,paste(aa),ii,"All"   ,] )*1000,
             type="l", lty=1, lwd=2,
             col=c("red","limegreen") )


###################################################
### code chunk number 19: Tr
###################################################
Tr <- list( "DM" = list( "Ins"       = DM.Ins,
                         "Dead"      = DM.Dead  ),
           "Ins" = list( "Dead(Ins)" = Ins.Dead ) )


###################################################
### code chunk number 20: make-ini
###################################################
str( ini <- Si[NULL,1:9] )


###################################################
### code chunk number 21: ini-fill
###################################################
ini[1:2,"lex.id"] <- 1:2
ini[1:2,"lex.Cst"] <- "DM"
ini[1:2,"Per"] <- 1995
ini[1:2,"Age"] <- 60
ini[1:2,"DMdur"] <- 5
ini[1:2,"sex"] <- c("M","F")
ini


###################################################
### code chunk number 22: simL
###################################################
set.seed(52381764)
Nsim <- 500
system.time( simL <- simLexis( Tr,
                              ini,
                          t.range = 12,
                                N = Nsim ) )


###################################################
### code chunk number 23: sum-simL
###################################################
summary( simL, by="sex" )


###################################################
### code chunk number 24: Tr.p-simP
###################################################
Tr.p <- list( "DM" = list( "Ins"       = DM.Ins,
                           "Dead"      = All.Dead  ),
             "Ins" = list( "Dead(Ins)" = All.Dead ) )
system.time( simP <- simLexis( Tr.p,
                                ini,
                            t.range = 12,
                                  N = Nsim ) )
summary( simP, by="sex" )


###################################################
### code chunk number 25: Cox-dur
###################################################
library( survival )
Cox.Dead <- coxph( Surv( DMdur, DMdur+lex.dur,
                         lex.Xst %in% c("Dead(Ins)","Dead")) ~
                   Ns( Age-DMdur, knots=ad.kn ) +
                   I(lex.Cst=="Ins") +
                   I(Per-2000) + sex,
               data = Si )
round( ci.exp( Cox.Dead ), 3 )


###################################################
### code chunk number 26: TR.c
###################################################
Tr.c <- list( "DM" = list( "Ins"       = Tr$DM$Ins,
                           "Dead"      = Cox.Dead  ),
             "Ins" = list( "Dead(Ins)" = Cox.Dead ) )
system.time( simC <- simLexis( Tr.c,
                                ini,
                            t.range = 12,
                                  N = Nsim ) )
summary( simC, by="sex" )


###################################################
### code chunk number 27: nState
###################################################
system.time(
nSt <- nState( subset(simL,sex=="M"),
               at=seq(0,11,0.2), from=1995, time.scale="Per" ) )
nSt[1:10,]


###################################################
### code chunk number 28: pstate0
###################################################
pM <- pState( nSt, perm=c(1,2,4,3) )
head( pM )
par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
plot( pM )
plot( pM, border="black", col="transparent", lwd=3 )
text( rep(as.numeric(rownames(pM)[nrow(pM)-1]),ncol(pM)),
      pM[nrow(pM),]-diff(c(0,pM[nrow(pM),]))/5,
      colnames( pM ), adj=1 )
box( col="white", lwd=3 )
box()


###################################################
### code chunk number 29: pstatex
###################################################
clr <- c("limegreen","orange")
# expand with a lighter version of the two chosen colors
clx <- c( clr, rgb( t( col2rgb( clr[2:1] )*2 + rep(255,3) ) / 3, max=255 ) )
par( mfrow=c(1,2), las=1, mar=c(3,3,4,2), mgp=c(3,1,0)/1.6 )
# Men
plot( pM, col=clx, xlab="Date of FU" )
lines( as.numeric(rownames(pM)), pM[,2], lwd=3 )
mtext( "60 year old male, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[2] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[1] )
axis( side=4 )
axis( side=4, at=1:19/20, labels=FALSE )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )
# Women
pF <- pState( nState( subset(simL,sex=="F"),
                      at=seq(0,11,0.2),
                      from=1995,
                      time.scale="Per" ),
              perm=c(1,2,4,3) )
plot( pF, col=clx, xlab="Date of FU" )
lines( as.numeric(rownames(pF)), pF[,2], lwd=3 )
mtext( "60 year old female, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[2] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[1] )
axis( side=4 )
axis( side=4, at=1:19/20, labels=FALSE )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )


###################################################
### code chunk number 30: pstatey
###################################################
par( mfrow=c(1,2), las=1, mar=c(3,3,4,2), mgp=c(3,1,0)/1.6 )
# Men
pM <- pState( nState( subset(simL,sex=="M"),
                      at=seq(0,11,0.2),
                      from=60,
                      time.scale="Age" ),
              perm=c(1,2,4,3) )
plot( pM, col=clx, xlab="Age" )
lines( as.numeric(rownames(pM)), pM[,2], lwd=3 )
mtext( "60 year old male, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[2] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[1] )
axis( side=4 )
axis( side=4, at=1:19/20, labels=FALSE )
axis( side=4, at=1:19/20, labels=FALSE, tcl=-0.4 )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )
# Women
pF <- pState( nState( subset(simL,sex=="F"),
                      at=seq(0,11,0.2),
                      from=60,
                      time.scale="Age" ),
              perm=c(1,2,4,3) )
plot( pF, col=clx, xlab="Age" )
lines( as.numeric(rownames(pF)), pF[,2], lwd=3 )
mtext( "60 year old female, diagnosed 1990, aged 55", side=3, line=2.5, adj=0, col=gray(0.6) )
mtext( "Survival curve", side=3, line=1.5, adj=0 )
mtext( "DM, no insulin   DM, Insulin", side=3, line=0.5, adj=0, col=clr[2] )
mtext( "DM, no insulin", side=3, line=0.5, adj=0, col=clr[1] )
axis( side=4 )
axis( side=4, at=1:9/10, labels=FALSE )
axis( side=4, at=1:19/20, labels=FALSE, tcl=-0.4 )
axis( side=4, at=1:99/100, labels=FALSE, tcl=-0.3 )


###################################################
### code chunk number 31: comp-0
###################################################
PrM  <- pState( nState( subset(simP,sex=="M"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
PrF  <- pState( nState( subset(simP,sex=="F"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
CoxM <- pState( nState( subset(simC,sex=="M"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
CoxF <- pState( nState( subset(simC,sex=="F"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )

par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(3,1,0)/1.6 )
 plot(   pM, border="black", col="transparent", lwd=3 )
lines(  PrM, border="blue" , col="transparent", lwd=3 )
lines( CoxM, border="red"  , col="transparent", lwd=3 )
text( 60.5, 0.05, "M" )
box( lwd=5, col="white" ) ; box( lwd=2, col="black" )

 plot(   pF, border="black", col="transparent", lwd=3 )
lines(  PrF, border="blue" , col="transparent", lwd=3 )
lines( CoxF, border="red"  , col="transparent", lwd=3 )
text( 60.5, 0.05, "F" )
box( lwd=5, col="white" ) ; box( lwd=2, col="black" )


###################################################
### code chunk number 32: 04simLexis.rnw:958-962
###################################################
ende <- Sys.time()
cat("  Start time:", format(anfang, "%F, %T"),
  "\n    End time:", format(  ende, "%F, %T"),
  "\nElapsed time:", round(difftime(ende, anfang, units = "mins"), 2), "minutes\n")


