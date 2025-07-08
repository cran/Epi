### R code from vignette source '03crisk.rnw'

###################################################
### code chunk number 1: 03crisk.rnw:28-40
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
### code chunk number 2: 03crisk.rnw:43-45
###################################################
anfang <- Sys.time()
cat("Start time:", format(anfang, "%F, %T"), "\n")


###################################################
### code chunk number 3: 03crisk.rnw:47-53
###################################################
vers <-
data.frame(R = substr(R.version.string, 11, 15),
         Epi = as.character(packageVersion(   "Epi")),
      popEpi = as.character(packageVersion("popEpi")))
names(vers) <- paste(" ", names(vers))
print(vers, row.names = FALSE)


###################################################
### code chunk number 4: 03crisk.rnw:251-266
###################################################
data(DMlate)
Ldm <- Lexis(entry = list(per = dodm,
                          age = dodm-dobth,
                          tfd = 0),
              exit = list(per = dox),
       exit.status = factor(!is.na(dodth), labels = c("DM", "Dead")),
              data = DMlate)
summary(Ldm, t = T)
set.seed(1952)
Mdm <- mcutLexis(Ldm,
                  wh = c('dooad','doins'),
          new.states = c('OAD','Ins'),
          seq.states = FALSE,
                ties = TRUE)
summary(Mdm)


###################################################
### code chunk number 5: 03crisk.rnw:271-276
###################################################
Sdm <- splitLexis(factorize(subset(Mdm,
                                   lex.Cst == "DM")),
                  time.scale = "tfd",
                      breaks = seq(0, 20, 1/12))
summary(Sdm)


###################################################
### code chunk number 6: boxes5
###################################################
boxes(Mdm, boxpos = list(x = c(15, 50, 15, 85, 85),
                         y = c(85, 50, 15, 85, 15)),
          scale.R = 100,
          show.BE = TRUE)


###################################################
### code chunk number 7: boxes4
###################################################
boxes(Relevel(Sdm, c(1, 4, 2, 3)),
      boxpos  = list(x = c(15, 85, 75, 15),
                     y = c(85, 85, 30, 15)),
      scale.R = 100,
      show.BE = TRUE )


###################################################
### code chunk number 8: 03crisk.rnw:314-317
###################################################
mD <- gamLexis(Sdm, ~ s(tfd, k = 5), to = 'Dead')
mO <- gamLexis(Sdm, ~ s(tfd, k = 5), to = 'OAD' )
mI <- gamLexis(Sdm, ~ s(tfd, k = 5), to = 'Ins' )


###################################################
### code chunk number 9: 03crisk.rnw:331-334
###################################################
nd <- data.frame(tfd = seq(0, 10, 1/20))
rownames(nd) <- nd$tfd
str(nd)


###################################################
### code chunk number 10: rates
###################################################
matshade(nd$tfd, cbind(ci.pred(mD, nd),
                       ci.pred(mI, nd),
                       ci.pred(mO, nd)) * 1000,
         col = c("black", "red", "blue"), log = "y", lwd = 3, plot = TRUE,
         xlab = "Time since DM diagnosis (years)",
         ylab = "Rates per 1000 PY", ylim = c(0.05, 500), yaxt = "n")
axis(side = 2, at = ll <- outer(c(1,2,5), -2:3, function(x,y) x * 10^y),
               labels = formatC(ll, digits = 4), las = 1)
axis(side = 2, at = outer(c(1.5,2:9), -2:3, function(x,y) x * 10^y),
               labels = NA, tcl = -0.3)
text(0, 0.5*0.6^c(1,2,0),
     c("Dead","Ins","OAD"),
     col = c("black","red","blue"), adj = 0, font = 2)


###################################################
### code chunk number 11: rates-l
###################################################
matshade(nd$tfd, cbind(ci.pred(mD, nd),
                       ci.pred(mI, nd),
                       ci.pred(mO, nd)) * 1000,
         col = c("black", "red", "blue"), lwd = 3, plot = TRUE,
         xlab = "Time since DM diagnosis (years)",
         ylab = "Rates per 1000 PY", ylim = c(0, 500), yaxs = "i")
text(8, 500 - c(2, 3, 1) * 20,
     c("Dead","Ins","OAD"),
     col = c("black","red","blue"), adj = 0, font = 2)


###################################################
### code chunk number 12: 03crisk.rnw:394-417
###################################################
# utility function to compute midpoints between sucessive values in a vector
mp <- function(x) x[-1] - diff(x) / 2
#
int <- 1 / 20
# rates at midpoints of intervals
lD <- mp(ci.pred(mD, nd)[, 1])
lI <- mp(ci.pred(mI, nd)[, 1])
lO <- mp(ci.pred(mO, nd)[, 1])
#
# cumulative rates and survival function at right border of the intervals
LD <- cumsum(lD) * int
LI <- cumsum(lI) * int
LO <- cumsum(lO) * int
# survival function, formula (1.1)
Sv <- exp(- LD - LI - LO)
#
# when integrating to get the cumulative *risks* we use the average
# of the survival function at the two endpoints
# (adding 1 as the first), formula (1.2)
Sv <- c(1, Sv)
rD <- c(0, cumsum(lD * mp(Sv)) * int)
rI <- c(0, cumsum(lI * mp(Sv)) * int)
rO <- c(0, cumsum(lO * mp(Sv)) * int)


###################################################
### code chunk number 13: 03crisk.rnw:422-426
###################################################
summary(rD + rI + rO + Sv)
oo <- options(digits = 20)
cbind(summary(Sv + rD + rI + rO))
options(oo)


###################################################
### code chunk number 14: stack
###################################################
zz <- mat2pol(cbind(rD, rI, rO, Sv), x = nd$tfd, # $
              xlim = c(0,10), xaxs = "i", yaxs = "i", las = 1,
              xlab = "Time since DM diagnosis (years)",
              ylab = "Probability",
               col =  c("black","red","blue","forestgreen"))
text(9, mp(zz["9", ]), c("Dead", "Ins", "OAD"," DM"), col = "white")
box(col = "white", lwd = 3)


###################################################
### code chunk number 15: 03crisk.rnw:460-465
###################################################
Sj <- c(sjA = sum(Sv * int),
        sjD = sum(rD * int),
        sjI = sum(rI * int),
        sjO = sum(rO * int))
c(Sj, sum(Sj))


###################################################
### code chunk number 16: 03crisk.rnw:515-517
###################################################
head(cbind(ci.pred(mI, nd),
           ci.exp (mI, nd)))


###################################################
### code chunk number 17: 03crisk.rnw:523-525
###################################################
str(ci.lin(mI, nd, sample = 4))
head(cbind(ci.pred(mI, nd), exp(ci.lin(mI, nd, sample = 4))))


###################################################
### code chunk number 18: 03crisk.rnw:571-578
###################################################
res <- ci.Crisk(list(OAD = mO,
                     Ins = mI,
                    Dead = mD),
                            nd = data.frame(tfd = seq(0, 10, 1/20)),
                            nB = 500,
                          perm = 4:1)
str(res)


###################################################
### code chunk number 19: 03crisk.rnw:613-620
###################################################
rsm <- ci.Crisk(list(OAD = mO,
                     Ins = mI,
                    Dead = mD),
                            nd = data.frame(tfd = seq(0, 10, 1/20)),
                            nB = 500,
                       sim.res = 'rates')
str(rsm)


###################################################
### code chunk number 20: 03crisk.rnw:628-635
###################################################
csm <- ci.Crisk(list(OAD = mO,
                     Ins = mI,
                    Dead = mD),
                            nd = data.frame(tfd = seq(0, 10, 1/20)),
                            nB = 500,
                       sim.res = 'crisk')
str(csm)


###################################################
### code chunk number 21: 03crisk.rnw:648-654
###################################################
Brates <- aperm(apply(rsm,
                      1:2,
                      quantile,
                      probs = c(.5, .025, .975)),
                c(2, 3, 1))
str(Brates)


###################################################
### code chunk number 22: rates-ci
###################################################
matshade(nd$tfd, cbind(ci.pred(mD, nd),
                       ci.pred(mI, nd),
                       ci.pred(mO, nd)) * 1000,
         ylim = c(0.1,500), yaxt = "n",
         ylab = "Rates per 1000 PY",
         xlab = "Time since DM diagnosis (years)",
         col = c("black", "red", "blue"), log = "y", lwd = 3, plot = TRUE)
matlines(nd$tfd,
         cbind(Brates[, "Dead", ],
               Brates[, "Ins" , ],
               Brates[, "OAD" , ]) * 1000,
         col = c("white", "black", "black"), lty = 3, lwd = c(3,1,1))
axis(side = 2, at = ll <- outer(c(1,2,5), -2:3, function(x, y) x * 10^y),
               labels = formatC(ll, digits = 4), las = 1)
axis(side = 2, at = outer(c(1.5, 2:9), -2:3, function(x, y) x * 10^y),
               labels = NA, tcl = -0.3)
text(0, 0.5 * 0.6^c(1,2,0),
     c("Dead", "Ins", "OAD"),
     col = c("black", "red", "blue"), adj = 0, font = 2)


###################################################
### code chunk number 23: crates
###################################################
matshade(res$time,
         cbind(res$Crisk[,"Dead",],
               res$Crisk[,"Ins" ,],
               res$Crisk[,"OAD" ,]), plot = TRUE,
         xlim = c(0,10), xaxs = "i", yaxs = "i", las = 1,
         xlab = "Time since DM diagnosis (years)",
         ylab = "Cumulative probability",
          col = c("black","red","blue"))
text(8, 0.3 + c(1, 0, 2) / 25,
     c("Dead", "Ins", "OAD"),
     col = c("black", "red", "blue"), adj = 0)


###################################################
### code chunk number 24: 03crisk.rnw:717-719
###################################################
str(res$Crisk)
str(res$Srisk)


###################################################
### code chunk number 25: stack-ci
###################################################
zz <- mat2pol(res$Crisk[,c("Dead", "Ins", "OAD", "Surv"),1],
              x = res$time,
           xlim = c(0, 10), xaxs = "i", yaxs = "i", las = 1,
           xlab = "Time since DM diagnosis (years)",
           ylab = "Probability",
            col =  c("black","red","blue","forestgreen") )
text(9, mp(zz["9",]), c("Dead", "Ins", "OAD", "DM"), col = "white" )
matshade(res$time,
         cbind(res$Srisk[, 1, ],
               res$Srisk[, 2, ],
               res$Srisk[, 3, ]),
         col = 'transparent', col.shade = "white", alpha = 0.4)


###################################################
### code chunk number 26: 03crisk.rnw:759-762
###################################################
s510 <- res$Stime[c("5", "10"),,]
dimnames(s510)[[1]] <- c(" 5 yr","10 yr")
round(ftable(s510, row.vars=1:2), 2)


###################################################
### code chunk number 27: 03crisk.rnw:778-781
###################################################
data(DMlate)
set.seed(7465)
wh <- sample(1:3, nrow(DMlate), replace = T, prob = c(4, 2, 6))


###################################################
### code chunk number 28: 03crisk.rnw:784-785
###################################################
wh[is.na(DMlate$dodth)] <- 0


###################################################
### code chunk number 29: 03crisk.rnw:790-792
###################################################
DMlate$codth <- factor(wh, labels = c("Alive", "CVD", "Can", "Oth"))
with(DMlate, table(codth, isDead = !is.na(dodth)))


###################################################
### code chunk number 30: 03crisk.rnw:801-803
###################################################
str(DMlate)
head(DMlate, 12)


###################################################
### code chunk number 31: 03crisk.rnw:810-817
###################################################
dmL <- Lexis(entry = list(per = dodm,
                          age = dodm - dobth,
                          tfD = 0),
              exit = list(per = dox),
       exit.status = codth,
              data = DMlate)
summary(dmL, t = T)


###################################################
### code chunk number 32: boxes
###################################################
boxes(dmL, boxpos = TRUE)


###################################################
### code chunk number 33: 03crisk.rnw:832-834
###################################################
sL <- splitLexis(dmL, time.scale = "age", breaks = seq(0, 120, 1/2))
summary(sL)


###################################################
### code chunk number 34: 03crisk.rnw:836-839
###################################################
mCVD <- gamLexis(sL, ~ s(tfD, by=sex), to = "CVD")
mCan <- gamLexis(sL, ~ s(tfD, by=sex), to = "Can")
mOth <- gamLexis(sL, ~ s(tfD, by=sex), to = "Oth")


###################################################
### code chunk number 35: 03crisk.rnw:841-844 (eval = FALSE)
###################################################
## mCVD <- glmLexis(sL, ~ Ns(tfD, kn=1:6*2):sex, to = "CVD")
## mCa  <- glmLexis(sL, ~ Ns(tfD, kn=1:6*2):sex, to = "Ca")
## mOth <- glmLexis(sL, ~ Ns(tfD, kn=1:6*2):sex, to = "Oth")


###################################################
### code chunk number 36: 03crisk.rnw:853-854
###################################################
nm <- data.frame(tfD = seq(0, 15, 1/20), sex = "M")


###################################################
### code chunk number 37: 03crisk.rnw:858-864
###################################################
cR <- ci.Crisk(list(CVD = mCVD,
                    Can = mCan,
                  Other = mOth),
               nB = 500,
               nd = nm)
str(cR)


###################################################
### code chunk number 38: cR
###################################################
clr <- c("black", "orange", "limegreen")
matshade(cR$time, cbind(cR$Crisk[, "CVD"  , ],
                        cR$Crisk[, "Can"  , ],
                        cR$Crisk[, "Other", ]),
         col = clr, lty = 1, lwd = 2,
         plot = TRUE, ylim = c(0, 1/3), yaxs = "i")
text(0, 1/3 - c(2,3,1)/30, c("CVD", "Can", "Oth"),
     col = clr, adj = 0, font = 2)


###################################################
### code chunk number 39: Sr1
###################################################
matshade(cR$time, cbind(cR$Srisk[,1,],
                        cR$Srisk[,2,],
                        cR$Srisk[,3,]),
         col = "black", lty = 1, lwd = 2,
         plot = TRUE, ylim = c(0,1), xaxs = "i", yaxs = "i")
text(14, mp(c(0, cR$Srisk["14", , 1], 1)),
     rev(c(dimnames(cR$Crisk)[[2]])))
box(bty = "o")


###################################################
### code chunk number 40: Sr2
###################################################
zz <- mat2pol(cR$Crisk[, c("Other", "Can", "CVD", "Surv"), "50%"],
              x = cR$time,
           xlim = c(0,15), xaxs = "i", yaxs = "i", las = 1,
           xlab = "Time since DM diagnosis (years)",
           ylab = "Probability",
            col =  c("gray", "red", "blue", "limegreen") )
matshade(cR$time, cbind(cR$Srisk[,1,],
                        cR$Srisk[,2,],
                        cR$Srisk[,3,]),
         col = "transparent", col.shade = "white", alpha = 0.4)
text(14, mp(c(0, cR$Srisk["14", , 1], 1)),
     rev(c(dimnames(cR$Crisk)[[2]])), col = "white")


###################################################
### code chunk number 41: 03crisk.rnw:949-950
###################################################
ftable(round(cR$Stime[paste(1:5 * 3), , ], 1), row.vars = 1)


###################################################
### code chunk number 42: 03crisk.rnw:973-992
###################################################
nm <- data.frame(tfD = seq(0, 15, 1/20), sex = "M")
nw <- data.frame(tfD = seq(0, 15, 1/20), sex = "F")
# set the seed
set.seed(1952)
mR <- ci.Crisk(list(CVD = mCVD,
                    Can = mCan,
                  Other = mOth),
               nd = nm,
               nB = 500,
          sim.res = "crisk" )
# REset the seed
set.seed(1952)
wR <- ci.Crisk(list(CVD = mCVD,
                    Can = mCan,
                  Other = mOth),
               nd = nw,
               nB = 500,
          sim.res = "crisk" )
str(wR)


###################################################
### code chunk number 43: 03crisk.rnw:997-1002
###################################################
dS <- mR[,"Surv",] - wR[,"Surv",]
dS <- apply(dS, 1, quantile, probs = c(.5, .025, .975)) * 100
str(dS)
rS <- mR[,"Surv",] / wR[,"Surv",]
rS <- apply(rS, 1, quantile, probs = c(.5, .025, .975))


###################################################
### code chunk number 44: difrat
###################################################
par(mfrow = c(1,2))
matshade(as.numeric(colnames(dS)), t(dS), plot = TRUE,
         lwd = 3, ylim = c(-5, 5),
         xlab = "Time since DM diagnosis (years)",
         ylab = "Men - Women survival difference (%)")
abline(h = 0)
matshade(as.numeric(colnames(rS)), t(rS), plot = TRUE,
         lwd = 3, ylim = c(1/1.2, 1.2), log ="y",
         xlab = "Time since DM diagnosis (years)",
         ylab = "Men - Women survival ratio")
abline(h = 1)


###################################################
### code chunk number 45: 03crisk.rnw:1028-1038
###################################################
fR <- ci.Crisk(list(CVD = mCVD,
                    Can = mCan,
                  Other = mOth),
               nd = nw,
               nB = 500,
          sim.res = "crisk" )
dxS <- mR[,"Surv",] - fR[,"Surv",]
dxS <- apply(dxS, 1, quantile, probs = c(.5, .025, .975)) * 100
rxS <- mR[,"Surv",] / fR[,"Surv",]
rxS <- apply(rxS, 1, quantile, probs = c(.5, .025, .975))


###################################################
### code chunk number 46: difratx
###################################################
par(mfrow = c(1,2))
matshade(as.numeric(colnames(dS)), t(dS), plot = TRUE,
         lwd = 3, ylim = c(-5, 5),
         xlab = "Time since DM diagnosis (years)",
         ylab = "Men - Women survival difference (%)")
matshade(as.numeric(colnames(dxS)), t(dxS), lty = 3, col = "forestgreen")
abline(h = 0)

matshade(as.numeric(colnames(rS)), t(rS), plot = TRUE,
         lwd = 3, ylim = c(1/1.2, 1.2), log ="y",
         xlab = "Time since DM diagnosis (years)",
         ylab = "Men - Women survival ratio")
matshade(as.numeric(colnames(rxS)), t(rxS), lty = 3, col = "forestgreen")
abline(h = 1)


###################################################
### code chunk number 47: 03crisk.rnw:1061-1065
###################################################
ende <- Sys.time()
cat("  Start time:", format(anfang, "%F, %T"),
  "\n    End time:", format(  ende, "%F, %T"),
  "\nElapsed time:", round(difftime(ende, anfang, units = "mins"), 2), "minutes\n")


