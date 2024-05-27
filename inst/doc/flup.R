### R code from vignette source 'flup.rnw'

###################################################
### code chunk number 1: flup.rnw:22-30
###################################################
options(width = 90, 
        SweaveHooks=list(fig=function()
                         par(mar = c(3, 3, 1, 1), 
                             mgp = c(3, 1, 0) / 1.6, 
                             las = 1,
                            lend = "butt",
                             bty = "n")))
library(Epi)


###################################################
### code chunk number 2: flup.rnw:33-35
###################################################
anfang <- Sys.time()
cat("Start time:", format(anfang, "%F, %T"), "\n")


###################################################
### code chunk number 3: flup.rnw:145-147
###################################################
library(Epi)
print( sessionInfo(), l = F)


###################################################
### code chunk number 4: flup.rnw:156-166
###################################################
data(DMlate)
head(DMlate)
dmL <- Lexis(entry = list(per = dodm, 
                          age = dodm-dobth, 
                          tfD = 0), 
              exit = list(per = dox), 
       exit.status = factor(!is.na(dodth), 
                            labels = c("DM", "Dead")), 
              data = DMlate)
timeScales(dmL)


###################################################
### code chunk number 5: flup.rnw:189-191
###################################################
str(dmL)
head(dmL)[, 1:10]


###################################################
### code chunk number 6: flup.rnw:207-208
###################################################
summary(dmL, timeScales = TRUE)


###################################################
### code chunk number 7: dmL1
###################################################
plot(dmL)


###################################################
### code chunk number 8: dmL2
###################################################
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0) / 1.6)
plot(dmL, 1:2, lwd = 1, col = c("blue", "red")[dmL$sex], 
     grid = TRUE, lty.grid = 1, col.grid = gray(0.7), 
     xlim = 1960 + c(0, 60), xaxs = "i", 
     ylim =   40 + c(0, 60), yaxs = "i", las = 1)
points(dmL, 1:2, pch = c(NA, 3)[dmL$lex.Xst], 
       col = "lightgray", lwd = 3, cex = 0.3)
points(dmL, 1:2, pch = c(NA, 3)[dmL$lex.Xst], 
       col = c("blue", "red")[dmL$sex], lwd = 1, cex = 0.3)
box(bty = 'o')


###################################################
### code chunk number 9: flup.rnw:264-267
###################################################
dmS1 <- splitLexis(dmL, "age", breaks = seq(0, 100, 5))
summary(dmL)
summary(dmS1)


###################################################
### code chunk number 10: flup.rnw:277-280
###################################################
wh.id <- c(9, 27, 52, 484)
subset(dmL , lex.id %in% wh.id)[, 1:10]
subset(dmS1, lex.id %in% wh.id)[, 1:10]


###################################################
### code chunk number 11: flup.rnw:286-288
###################################################
dmS2 <- splitLexis(dmS1, "tfD", breaks = c(0, 1, 2, 5, 10, 20, 30, 40))
subset(dmS2, lex.id %in% wh.id)[, 1:10]


###################################################
### code chunk number 12: flup.rnw:293-303
###################################################
if (require(popEpi, quietly = TRUE)) 
   {
   options("popEpi.datatable" = FALSE)
   dmM <- splitMulti(dmL, 
                     age = seq(0, 100, 5), 
                     tfD = c(0, 1, 2, 5, 10, 20, 30, 40), 
                    drop = FALSE)
   summary(dmS2)
   summary(dmM)
   }


###################################################
### code chunk number 13: flup.rnw:314-320
###################################################
if (require(popEpi, quietly = TRUE)) 
   {
   identical(dmS2, dmM)
   class(dmS2)
   class(dmM)
   }


###################################################
### code chunk number 14: flup.rnw:350-357
###################################################
subset(dmL, lex.id %in% wh.id)
dmC <- cutLexis(data = dmL, 
                 cut = dmL$doins, 
           timescale = "per", 
           new.state = "Ins", 
           new.scale = "tfI")
subset(dmC, lex.id %in% wh.id)[, 1:10]


###################################################
### code chunk number 15: flup.rnw:372-378
###################################################
dmS2C <- cutLexis(data = dmS2, 
                   cut = dmS2$doins, 
             timescale = "per", 
             new.state = "Ins", 
             new.scale = "tfI")
subset(dmS2C, lex.id %in% wh.id)


###################################################
### code chunk number 16: flup.rnw:403-404
###################################################
summary(dmS2C, timeScales = TRUE)


###################################################
### code chunk number 17: box1
###################################################
boxes(dmC, boxpos = TRUE, scale.R = 1000, show.BE = TRUE)


###################################################
### code chunk number 18: flup.rnw:447-455
###################################################
timeBand(dmS2C, "age", "middle")[1:10]
# For nice printing and column labelling we use the data.frame() function:
data.frame(dmS2C[, c("per", "age", "tfD", "lex.dur")], 
           mid.age = timeBand(dmS2C, "age", "middle"), 
             mid.t = timeBand(dmS2C, "tfD", "middle"), 
            left.t = timeBand(dmS2C, "tfD", "left"  ), 
           right.t = timeBand(dmS2C, "tfD", "right" ), 
            fact.t = timeBand(dmS2C, "tfD", "factor"))[1:15, ]


###################################################
### code chunk number 19: flup.rnw:491-492
###################################################
summary((dmS2$age - dmS2$tfD) - (dmS2$dodm - dmS2$dobth))


###################################################
### code chunk number 20: flup.rnw:497-500
###################################################
summary(timeBand(dmS2, "age", "middle") -
        timeBand(dmS2, "tfD", "middle") - 
        (dmS2$dodm - dmS2$dobth))


###################################################
### code chunk number 21: flup.rnw:605-607
###################################################
dmCs <- splitLexis(dmC, time.scale = "age", breaks = seq(0, 110, 1/4))
summary(dmCs, t = T)


###################################################
### code chunk number 22: flup.rnw:629-634
###################################################
(a.kn <- with(subset(dmCs, lex.Xst == "Dead"), 
              quantile(age+lex.dur, (1:5-0.5)/5)))
(i.kn <- c(0, 
           with(subset(dmCs, lex.Xst == "Dead" & lex.Cst == "Ins"), 
                quantile(tfI+lex.dur, (1:4)/5))))


###################################################
### code chunk number 23: flup.rnw:650-655
###################################################
ma <- glm((lex.Xst == "Dead") ~ Ns(age, knots = a.kn), 
           family = poisson, 
           offset = log(lex.dur), 
             data = dmCs)
summary(ma)


###################################################
### code chunk number 24: flup.rnw:674-678
###################################################
Ma <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(age, knots = a.kn), 
          family = poisreg, 
            data = dmCs)
summary(Ma)


###################################################
### code chunk number 25: flup.rnw:686-688
###################################################
Xa <- glm.Lexis(dmCs, from = "DM", to = "Dead", 
                formula = ~ Ns(age, knots = a.kn))


###################################################
### code chunk number 26: flup.rnw:691-692
###################################################
attr(Xa, "Lexis")


###################################################
### code chunk number 27: flup.rnw:701-702
###################################################
xa <- glm.Lexis(dmCs, formula = ~ Ns(age, knots = a.kn))


###################################################
### code chunk number 28: flup.rnw:705-706
###################################################
c(deviance(ma), deviance(Ma), deviance(Xa), deviance(xa))


###################################################
### code chunk number 29: pr-a
###################################################
nd <- data.frame(age = 40:85, lex.dur = 1000)
pr.0 <- ci.pred(ma, newdata = nd)      # mortality per 100 PY
pr.a <- ci.pred(Ma, newdata = nd)*1000 # mortality per 100 PY
summary(pr.0/pr.a)
matshade(nd$age, pr.a, plot = TRUE, 
          type = "l", lty = 1, 
          log = "y", xlab = "Age (years)", 
          ylab = "DM mortality per 1000 PY")


###################################################
### code chunk number 30: flup.rnw:753-758
###################################################
pm <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(age, knots = a.kn) 
                                              + lex.Cst + sex, 
          family = poisreg, 
            data = dmCs)
round(ci.exp(pm), 3)


###################################################
### code chunk number 31: flup.rnw:772-777
###################################################
pm <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(age, knots = a.kn) 
                                            + Ns(tfI, knots = i.kn) 
                                            + lex.Cst + sex, 
          family = poisreg, 
            data = tsNA20(dmCs))


###################################################
### code chunk number 32: flup.rnw:783-789
###################################################
Pm <- glm.Lexis(tsNA20(dmCs), 
                form = ~ Ns(age, knots = a.kn) 
                       + Ns(tfI, knots = i.kn) 
                       + lex.Cst + sex)
c(deviance(Pm), deviance(pm))
identical(model.matrix(Pm), model.matrix(pm))


###################################################
### code chunk number 33: flup.rnw:795-796
###################################################
round(ci.exp(Pm, subset = "ex"), 3)


###################################################
### code chunk number 34: ins-time
###################################################
ndI <- data.frame(expand.grid(tfI = c(NA, seq(0, 15, 0.1)), 
                               ai = seq(40, 80, 10)), 
                  sex = "M", 
                  lex.Cst = "Ins")
ndI <- transform(ndI, age = ai+tfI)
head(ndI)
ndA <- data.frame(age = seq(40, 100, 0.1), tfI = 0,  lex.Cst = "DM", sex = "M")
pri <- ci.pred(Pm, ndI) * 1000
pra <- ci.pred(Pm, ndA) * 1000
matshade(ndI$age, pri, plot = TRUE, las = 1, 
         xlab = "Age (years)", ylab = "DM mortality per 1000 PY", 
         log = "y", lty = 1, col = "blue")
matshade(ndA$age, pra)


###################################################
### code chunk number 35: flup.rnw:833-837
###################################################
library(survival)
cm <- coxph(Surv(age, age+lex.dur, lex.Xst == "Dead") ~
            Ns(tfI, knots = i.kn) + lex.Cst + sex, 
            data = tsNA20(dmCs))


###################################################
### code chunk number 36: flup.rnw:841-844
###################################################
Cm <- coxph.Lexis(tsNA20(dmCs), 
                  form = age ~ Ns(tfI, knots = i.kn) + lex.Cst + sex)
cbind(ci.exp(cm), ci.exp(Cm))


###################################################
### code chunk number 37: flup.rnw:853-856
###################################################
round(cbind(ci.exp(Pm), 
       rbind(matrix(NA, 5, 3), 
             ci.exp(cm)[-6, ])), 3)


###################################################
### code chunk number 38: Ieff
###################################################
nd <- data.frame(tfI = seq(0, 15, , 151), lex.Cst = "Ins", sex = "M")
nr <- data.frame(tfI =     2            , lex.Cst = "Ins", sex = "M")
ppr <- ci.exp(pm, list(nd, nr), xvars = "age")
cpr <- ci.exp(cm, list(nd, nr))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(nd$tfI, cbind(ppr, cpr), plot = T, 
         lty = c(1, 2), log = "y", 
         xlab = "Time since insulin (years)", ylab = "Rate ratio")
abline(h = 1, lty = 3)


###################################################
### code chunk number 39: IeffR
###################################################
nd <- data.frame(tfI = seq(0, 15, , 151), lex.Cst = "Ins", sex = "M")
nr <- data.frame(tfI =     0            , lex.Cst = "DM" , sex = "M")
ppr <- ci.exp(pm, list(nd, nr), xvars = "age")
cpr <- ci.exp(cm, list(nd, nr))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(nd$tfI, cbind(ppr, cpr), 
         xlab = "Time since insulin (years)", 
         ylab = "Rate ratio relative to non-Insulin", 
         lty = c(1, 2), log = "y", plot = T)


###################################################
### code chunk number 40: flup.rnw:963-968
###################################################
imx <- glm.Lexis(tsNA20(dmCs), 
                 formula = ~ Ns(age      , knots = a.kn) 
                           + Ns(      tfI, knots = i.kn)
                           + Ns(age - tfI, knots = a.kn)
                           + lex.Cst + sex)


###################################################
### code chunk number 41: flup.rnw:978-988
###################################################
Im <- glm.Lexis(tsNA20(dmCs), 
                formula = ~ Ns(age       , knots = a.kn) 
                          + Ns(       tfI, knots = i.kn)
                          + Ns((age - tfI) * (lex.Cst == "Ins"), knots = a.kn)
                          + lex.Cst + sex)
im <- glm.Lexis(tsNA20(dmCs), 
                formula = ~ Ns(age      , knots = a.kn) 
                          + Ns(      tfI, knots = i.kn)
                  + lex.Cst:Ns(age - tfI, knots = a.kn)
                          + lex.Cst + sex)


###################################################
### code chunk number 42: flup.rnw:1003-1004
###################################################
anova(imx, Im, im, test = 'Chisq')


###################################################
### code chunk number 43: dur-int
###################################################
pxi <- ci.pred(imx, ndI)
pxa <- ci.pred(imx, ndA)
pIi <- ci.pred(Im , ndI)
pIa <- ci.pred(Im , ndA)
pii <- ci.pred(im , ndI)
pia <- ci.pred(im , ndA)
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, cbind(pxi, pIi, pii)*1000, plot = T, log = "y", 
         xlab = "Age", ylab = "Mortality per 1000 PY", 
         lty = 1, lwd = 2, col = c("blue", "forestgreen", "red"), alpha = 0.1)
matshade(ndA$age, cbind(pxa, pIa, pia)*1000, 
         lty = 1, lwd = 2, col = c("blue", "forestgreen", "red"), alpha = 0.1)


###################################################
### code chunk number 44: dur-int-RR
###################################################
ndR <- transform(ndI, tfI = 0, lex.Cst = "DM")
cbind(head(ndI), head(ndR))
Rxi <- ci.exp(imx, list(ndI, ndR))
Rii <- ci.exp(im , list(ndI, ndR))
RIi <- ci.exp(Im , list(ndI, ndR))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, cbind(Rxi, RIi, Rii), plot = T, log = "y", 
         xlab = "Age (years)", ylab = "Rate ratio vs, non-Insulin", 
         lty = 1, lwd = 2, col = c("blue", "forestgreen", "red"), alpha = 0.1)
abline(h = 1)
abline(h = ci.exp(imx, subset = "lex.Cst")[, 1], lty = "25", col = "blue")


###################################################
### code chunk number 45: splint
###################################################
gm <- glm.Lexis(tsNA20(dmCs), 
                formula = ~ Ns(age, knots = a.kn) 
                          + Ns(tfI, knots = i.kn)
                          + lex.Cst:Ns(age, knots = a.kn):Ns(tfI, knots = i.kn)
                          + lex.Cst + sex)
pgi <- ci.pred(gm, ndI)
pga <- ci.pred(gm, ndA)
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, cbind(pgi, pii)*1000,  plot = T, 
         lty = c("solid", "21"), lend = "butt", lwd = 2, log = "y", 
         xlab = "Age (years)", ylab = "Mortality rates per 1000 PY", 
         alpha = c(0.2, 0.1), col = c("black", "red"))
matshade(ndA$age, cbind(pga, pia)*1000, 
         lty = c("solid", "21"), lend = "butt", lwd = 2, 
         alpha = c(0.2, 0.1), col = c("black", "red"))


###################################################
### code chunk number 46: RR-int
###################################################
ndR <- transform(ndI, lex.Cst = "DM", tfI = 0)
iRR <- ci.exp(im, ctr.mat = list(ndI, ndR))
gRR <- ci.exp(gm, ctr.mat = list(ndI, ndR))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, cbind(gRR, iRR), lty = 1, log = "y", plot = TRUE, 
         xlab = "Age (years)", ylab = "Rate ratio: Ins vs. non-Ins", 
         col = c("black", "red"))
abline(h = 1)


###################################################
### code chunk number 47: flup.rnw:1127-1140
###################################################
dmd <- glm.Lexis(dmCs, 
                 from = "DM", to = "Dead", 
                 formula = ~ Ns(age, knots = a.kn) 
                           + sex)
ind <- glm.Lexis(dmCs, 
                 from = "Ins", to = "Dead", 
                 formula = ~ Ns(age      , knots = a.kn) 
                           + Ns(      tfI, knots = i.kn)
                           + Ns(age - tfI, knots = a.kn)
                           + sex)
ini <- ci.pred(ind, ndI)
dmi <- ci.pred(dmd, ndI)
dma <- ci.pred(dmd, ndA)


###################################################
### code chunk number 48: sep-mort
###################################################
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, ini*1000, plot = TRUE, log = "y", 
         xlab = "Age (years)", ylab = "Mortality rates per 1000 PY", 
         lwd = 2, col = "red")
matshade(ndA$age, dma*1000, 
         lwd = 2, col = "black")


###################################################
### code chunk number 49: sep-HR
###################################################
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, ci.ratio(ini, dmi), plot = TRUE, log = "y", 
         xlab = "Age (years)", ylab = "RR insulin vs. no insulin", 
         lwd = 2, col = "red")
abline(h = 1)


###################################################
### code chunk number 50: flup.rnw:1186-1193
###################################################
dmCs <- cutLexis(data = dmS2, 
                  cut = dmS2$doins, 
            timescale = "per", 
            new.state = "Ins", 
            new.scale = "tfI", 
         split.states = TRUE)
summary(dmCs)


###################################################
### code chunk number 51: box4
###################################################
boxes(dmCs, boxpos = list(x = c(15, 15, 85, 85), 
                          y = c(85, 15, 85, 15)), 
      scale.R = 1000, show.BE = TRUE)


###################################################
### code chunk number 52: flup.rnw:1223-1230
###################################################
dmM <- mcutLexis(dmL, 
           timescale = "per", 
                  wh = c("doins", "dooad"), 
          new.states = c("Ins", "OAD"), 
          new.scales = c("tfI", "tfO"), 
        ties.resolve = TRUE)
summary(dmM, t = T)


###################################################
### code chunk number 53: flup.rnw:1234-1237
###################################################
wh <- c(subset(dmM, lex.Cst == "Ins-OAD")$lex.id[1:2], 
        subset(dmM, lex.Cst == "OAD-Ins")$lex.id[1:2])
subset(dmM, lex.id %in% wh)


###################################################
### code chunk number 54: mbox
###################################################
boxes(dmM, boxpos = list(x = c(15, 80, 40, 40, 85, 85), 
                         y = c(50, 50, 90, 10, 90, 10)), 
           scale.R = 1000, show.BE = TRUE)


###################################################
### code chunk number 55: mboxr
###################################################
summary(dmMr <- Relevel(dmM, list('OAD+Ins' = 5:6), first = FALSE))
boxes(dmMr, boxpos = list(x = c(15, 50, 15, 85, 85), 
                          y = c(85, 50, 15, 85, 15)), 
            scale.R = 1000, show.BE = TRUE)


###################################################
### code chunk number 56: flup.rnw:1406-1410
###################################################
ende <- Sys.time()
cat("  Start time:", format(anfang, "%F, %T"), "\n")
cat("    End time:", format(  ende, "%F, %T"), "\n")
cat("Elapsed time:", round(difftime(ende, anfang, units = "mins"), 2), "minutes\n")


