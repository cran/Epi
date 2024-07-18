### R code from vignette source 'yll.rnw'

###################################################
### code chunk number 1: yll.rnw:21-27
###################################################
options(width=90,
         SweaveHooks=list(fig=function()
         par(mar = c(3, 3, 1, 1),
             mgp = c(3, 1, 0) / 1.6,
             las = 1,
             bty = "n")))


###################################################
### code chunk number 2: yll.rnw:31-33
###################################################
anfang <- Sys.time()
cat("Start time:", format(anfang, "%F, %T"), "\n")


###################################################
### code chunk number 3: states
###################################################
library(Epi)
TM <- matrix(NA, 4, 4)
rownames(TM) <-
colnames(TM) <- c("Well", "DM", "Dead", "Dead(DM)")
TM[1, 2:3] <- TM[2, 4] <- 1
TM
zz <- boxes(TM, boxpos = list(x = c(20, 80, 20, 80),
                              y = c(80, 80, 20, 20)),
                wm = 1.5,
                hm = 4)


###################################################
### code chunk number 4: states
###################################################
zz$Arrowtext <- c(expression(lambda(a)),
                  expression(mu[W](a)),
                  expression(mu[D][M](a,d)))
boxes.MS(zz)


###################################################
### code chunk number 5: yll.rnw:280-281
###################################################
data(DMepi)


###################################################
### code chunk number 6: yll.rnw:287-289
###################################################
str(DMepi)
head(DMepi)


###################################################
### code chunk number 7: yll.rnw:309-315
###################################################
DMepi <- transform(subset(DMepi, A > 30),
                   A = A + 0.5,
                   P = P + 0.5,
                 D.T = D.nD + D.DM,
                 Y.T = Y.nD + Y.DM)
head(DMepi)


###################################################
### code chunk number 8: yll.rnw:321-346
###################################################
# Knots used in all models
(a.kn <- seq(40, 95, , 6))
(p.kn <- seq(1997, 2015, , 4))
(c.kn <- seq(1910, 1976, , 6))
# Check the number of events between knots
ae <- xtabs(cbind(D.nD, D.DM, X) ~ cut(A, c(30, a.kn, Inf)) + sex, data=DMepi)
ftable(addmargins(ae, 1), col.vars=3:2)
pe <- xtabs(cbind(D.nD, D.DM, X) ~ cut(P, c(1990, p.kn, Inf)) + sex, data=DMepi)
ftable(addmargins(pe, 1), col.vars=3:2)
ce <- xtabs(cbind(D.nD, D.DM, X) ~ cut(P-A, c(-Inf, c.kn, Inf)) + sex, data=DMepi)
ftable(addmargins(ce, 1), col.vars=3:2)
# Fit an APC-model for all transitions, separately for men and women
mW.m <- glm(cbind(D.nD, Y.nD) ~ -1 + Ns(    A, knots=a.kn, int=TRUE) +
                                     Ns(P    , knots=p.kn, ref=2005) +
                                     Ns(P - A, knots=c.kn, ref=1950),
            family = poisreg,
              data = subset(DMepi, sex=="M"))
mD.m <- update(mW.m, cbind(D.DM, Y.DM) ~ .)
mT.m <- update(mW.m, cbind(D.T , Y.T ) ~ .)
lW.m <- update(mW.m, cbind(X   , Y.nD) ~ .)
# Model for women
mW.f <- update(mW.m, data = subset(DMepi, sex == "F"))
mD.f <- update(mD.m, data = subset(DMepi, sex == "F"))
mT.f <- update(mT.m, data = subset(DMepi, sex == "F"))
lW.f <- update(lW.m, data = subset(DMepi, sex == "F"))


###################################################
### code chunk number 9: yll.rnw:353-390
###################################################
a.ref <- 30:90
p.ref <- 1996:2016
aYLL <- NArray(list(type = c("Imm", "Tot", "Sus"),
                       sex = levels(DMepi$sex),
                       age = a.ref,
                      date = p.ref))
str(aYLL)
system.time(
for(ip in p.ref)
   {
   nd <- data.frame(A = seq(30, 90, 0.2)+0.1,
                     P = ip,
                  Y.nD = 1,
                  Y.DM = 1,
                  Y.T  = 1)
   muW.m <- ci.pred(mW.m, nd)[, 1]
   muD.m <- ci.pred(mD.m, nd)[, 1]
   muT.m <- ci.pred(mT.m, nd)[, 1]
   lam.m <- ci.pred(lW.m, nd)[, 1]
   muW.f <- ci.pred(mW.f, nd)[, 1]
   muD.f <- ci.pred(mD.f, nd)[, 1]
   muT.f <- ci.pred(mT.f, nd)[, 1]
   lam.f <- ci.pred(lW.f, nd)[, 1]
   aYLL["Imm", "M", , paste(ip)] <- yll(int=0.2, muW.m, muD.m, lam=NULL,
                                      A=a.ref, age.in=30, note=FALSE)[-1]
   aYLL["Imm", "F", , paste(ip)] <- yll(int=0.2, muW.f, muD.f, lam=NULL,
                                      A=a.ref, age.in=30, note=FALSE)[-1]
   aYLL["Tot", "M", , paste(ip)] <- yll(int=0.2, muT.m, muD.m, lam=NULL,
                                      A=a.ref, age.in=30, note=FALSE)[-1]
   aYLL["Tot", "F", , paste(ip)] <- yll(int=0.2, muT.f, muD.f, lam=NULL,
                                      A=a.ref, age.in=30, note=FALSE)[-1]
   aYLL["Sus", "M", , paste(ip)] <- yll(int=0.2, muW.m, muD.m, lam=lam.m,
                                      A=a.ref, age.in=30, note=FALSE)[-1]
   aYLL["Sus", "F", , paste(ip)] <- yll(int=0.2, muW.f, muD.f, lam=lam.f,
                                      A=a.ref, age.in=30, note=FALSE)[-1]
   })
round(ftable(aYLL[, , seq(1, 61, 10), ], col.vars=c(3, 2)), 1)


###################################################
### code chunk number 10: imm
###################################################
plyll <- function(wh, xtxt){
par(mfrow = c(1, 2),
      mar = c(3, 3, 1, 1),
      mgp = c(3, 1, 0) / 1.6,
      bty = "n",
      las = 1)
matplot(a.ref, aYLL[wh, "M", , ],
         type="l", lty=1, col="blue", lwd=1:2,
         ylim=c(0, 12), xlab="Age",
         ylab=paste0("Years lost to DM", xtxt),
         yaxs="i")
abline(v=50, h=1:11, col=gray(0.7))
text(90, 11.5, "Men", col="blue", adj=1)
text(40, aYLL[wh, "M", "40", "1996"], "1996", adj=c(0, 0), col="blue")
text(43, aYLL[wh, "M", "44", "2016"], "2016", adj=c(1, 1), col="blue")

matplot(a.ref, aYLL[wh, "F", , ],
         type="l", lty=1, col="red", lwd=1:2,
         ylim=c(0, 12), xlab="Age",
         ylab=paste0("Years lost to DM", xtxt),
         yaxs="i")
abline(v=50, h=1:11, col=gray(0.7))
text(90, 11.5, "Women", col="red", adj=1)
text(40, aYLL[wh, "F", "40", "1996"], "1996", adj=c(0, 0), col="red")
text(43, aYLL[wh, "F", "44", "2016"], "2016", adj=c(1, 1), col="red")
}
plyll("Imm", " - immunity assumption")


###################################################
### code chunk number 11: tot
###################################################
plyll("Tot", " - total mortality refernce")


###################################################
### code chunk number 12: sus
###################################################
plyll("Sus", " - susceptibility assumed")


###################################################
### code chunk number 13: yll.rnw:473-477
###################################################
ende <- Sys.time()
cat("  Start time:", format(anfang, "%F, %T"), "\n")
cat("    End time:", format(  ende, "%F, %T"), "\n")
cat("Elapsed time:", round(difftime(ende, anfang, units = "mins"), 2), "minutes\n")


