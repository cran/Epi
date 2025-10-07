### R code from vignette source '01flup.rnw'

###################################################
### code chunk number 1: 01flup.rnw:29-41
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
### code chunk number 2: 01flup.rnw:44-46
###################################################
anfang <- Sys.time()
cat("Start time:", format(anfang, "%F, %T"), "\n")


###################################################
### code chunk number 3: 01flup.rnw:48-54
###################################################
vers <-
data.frame(R = substr(R.version.string, 11, 15),
         Epi = as.character(packageVersion(   "Epi")),
      popEpi = as.character(packageVersion("popEpi")))
names(vers) <- paste(" ", names(vers))
print(vers, row.names = FALSE)


###################################################
### code chunk number 4: 01flup.rnw:223-233
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
### code chunk number 5: 01flup.rnw:256-258
###################################################
str(dmL)
head(dmL)[, 1:11]


###################################################
### code chunk number 6: 01flup.rnw:275-276
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
### code chunk number 9: 01flup.rnw:334-337
###################################################
dmS1 <- splitLexis(dmL, "age", breaks = seq(0, 100, 5))
summary(dmL)
summary(dmS1)


###################################################
### code chunk number 10: 01flup.rnw:347-350
###################################################
wh.id <- c(9, 27, 52, 484)
subset(dmL , lex.id %in% wh.id)[, 1:10]
subset(dmS1, lex.id %in% wh.id)[, 1:10]


###################################################
### code chunk number 11: 01flup.rnw:358-360
###################################################
dmS2 <- splitLexis(dmS1, "tfD", breaks = c(0, 1, 2, 5, 10, 20, 30, 40))
subset(dmS2, lex.id %in% wh.id)[, 1:10]


###################################################
### code chunk number 12: 01flup.rnw:365-371
###################################################
dmM <- splitMulti(dmL,
                  age = seq(0, 100, 5),
                  tfD = c(0, 1, 2, 5, 10, 20, 30, 40),
                 drop = FALSE)
summary(dmS2)
summary(dmM)


###################################################
### code chunk number 13: 01flup.rnw:403-410
###################################################
subset(dmL, lex.id %in% wh.id)[, 1:11]
dmC <- cutLexis(data = dmL,
                 cut = dmL$doins,
           timescale = "per",
           new.state = "Ins",
           new.scale = "tfI")
subset(dmC, lex.id %in% wh.id)[, 1:11]


###################################################
### code chunk number 14: 01flup.rnw:425-431
###################################################
dmS2C <- cutLexis(data = dmS2,
                   cut = dmS2$doins,
             timescale = "per",
             new.state = "Ins",
             new.scale = "tfI")
subset(dmS2C, lex.id %in% wh.id)[, 1:11]


###################################################
### code chunk number 15: 01flup.rnw:449-450
###################################################
summary(dmS2C, timeScales = TRUE)


###################################################
### code chunk number 16: box1
###################################################
boxes(dmC, boxpos = TRUE, scale.R = 1000, show.BE = TRUE)
legendbox(70, 95)


###################################################
### code chunk number 17: 01flup.rnw:495-503
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
### code chunk number 18: 01flup.rnw:534-535
###################################################
summary((dmS2$age - dmS2$tfD) - (dmS2$dodm - dmS2$dobth))


###################################################
### code chunk number 19: 01flup.rnw:541-544
###################################################
summary(timeBand(dmS2, "age", "middle") -
        timeBand(dmS2, "tfD", "middle") -
        (dmS2$dodm - dmS2$dobth))


###################################################
### code chunk number 20: 01flup.rnw:653-655
###################################################
dmCs <- splitLexis(dmC, time.scale = "age", breaks = seq(0, 110, 1/4))
summary(dmCs, t = T)


###################################################
### code chunk number 21: 01flup.rnw:677-682
###################################################
(a.kn <- with(subset(dmCs, lex.Xst == "Dead"),
              quantile(age+lex.dur, seq(5, 95, , 5)  /100)))
(i.kn <- c(0,
           with(subset(dmCs, lex.Xst == "Dead" & lex.Cst == "Ins"),
                quantile(tfI+lex.dur, seq(20, 95, , 4) / 100))))


###################################################
### code chunk number 22: 01flup.rnw:698-703
###################################################
ma <- glm((lex.Xst == "Dead") ~ Ns(age, knots = a.kn),
           family = poisson,
           offset = log(lex.dur),
             data = dmCs)
summary(ma)


###################################################
### code chunk number 23: 01flup.rnw:722-726
###################################################
Ma <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(age, knots = a.kn),
          family = poisreg,
            data = dmCs)
summary(Ma)


###################################################
### code chunk number 24: 01flup.rnw:732-734
###################################################
Xa <- glmLexis(dmCs, formula = ~ Ns(age, knots = a.kn),
                     from = "DM", to = "Dead",)


###################################################
### code chunk number 25: 01flup.rnw:739-740
###################################################
attr(Xa, "Lexis")


###################################################
### code chunk number 26: 01flup.rnw:751-754
###################################################
transient(dmCs)
absorbing(dmCs)
preceding(dmCs, absorbing(dmCs))


###################################################
### code chunk number 27: 01flup.rnw:758-759
###################################################
xa <- glmLexis(dmCs, formula = ~ Ns(age, knots = a.kn))


###################################################
### code chunk number 28: 01flup.rnw:762-766
###################################################
c(ma = deviance(ma),
  Ma = deviance(Ma),
  Xa = deviance(Xa),
  xa = deviance(xa))


###################################################
### code chunk number 29: pr-a
###################################################
nd <- data.frame(age = 40:85, lex.dur = 1000)
pr.0 <- ci.pred(ma, newdata = nd)      # mortality per 1000 PY
pr.a <- ci.pred(Ma, newdata = nd)*1000 # mortality per 1000 PY
summary(pr.0 / pr.a)
matshade(nd$age, pr.a, plot = TRUE,
         type = "l", lty = 1,
         log = "y", xlab = "Age (years)",
         ylab = "DM mortality per 1000 PY")


###################################################
### code chunk number 30: 01flup.rnw:813-818
###################################################
pm <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(age, knots = a.kn)
                                              + lex.Cst + sex,
          family = poisreg,
            data = dmCs)
round(ci.exp(pm), 3)


###################################################
### code chunk number 31: 01flup.rnw:821-823
###################################################
pm <- glmLexis(dmCs, ~ Ns(age, knots = a.kn) + lex.Cst + sex)
round(ci.exp(pm), 3)


###################################################
### code chunk number 32: 01flup.rnw:840-845
###################################################
pm <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(age, knots = a.kn)
                                            + Ns(tfI, knots = i.kn)
                                            + lex.Cst + sex,
          family = poisreg,
            data = tsNA20(dmCs))


###################################################
### code chunk number 33: 01flup.rnw:851-857
###################################################
Pm <- glmLexis(tsNA20(dmCs),
               form = ~ Ns(age, knots = a.kn)
                      + Ns(tfI, knots = i.kn)
                      + lex.Cst + sex)
c(deviance(Pm), deviance(pm))
identical(model.matrix(Pm), model.matrix(pm))


###################################################
### code chunk number 34: 01flup.rnw:863-864
###################################################
round(ci.exp(Pm, subset = "ex"), 3)


###################################################
### code chunk number 35: ins-time
###################################################
ndI <- data.frame(expand.grid(tfI = c(NA, seq(0, 15, 0.1)),
                               ai = seq(40, 80, 10)),
                  lex.Cst = "Ins",
                      sex = "M")
ndI <- transform(ndI, age = ai + tfI)
head(ndI)
ndA <- data.frame(age = seq(40, 100, 0.1),
                  tfI = 0,
              lex.Cst = "DM",
                  sex = "M")
pri <- ci.pred(Pm, ndI) * 100
pra <- ci.pred(Pm, ndA) * 100
matshade(ndI$age, pri, plot = TRUE,
         xlab = "Attained age (years)", ylab = "DM mortality per 100 PY",
         las = 1, log = "y", lty = 1, col = "blue")
matshade(ndA$age, pra)


###################################################
### code chunk number 36: 01flup.rnw:933-936
###################################################
cm <- coxph(Surv(age, age + lex.dur, lex.Xst == "Dead") ~
            Ns(tfI, knots = i.kn) + lex.Cst + sex,
            data = tsNA20(dmCs))


###################################################
### code chunk number 37: 01flup.rnw:940-943
###################################################
Cm <- coxphLexis(tsNA20(dmCs),
                  formula = age ~ Ns(tfI, knots = i.kn) + lex.Cst + sex)
round(cbind(ci.exp(cm), ci.exp(Cm)), 4)


###################################################
### code chunk number 38: 01flup.rnw:960-963
###################################################
round(cbind(ci.exp(Pm),
       rbind(matrix(NA, 5, 3),
             ci.exp(cm)[-6, ])), 3)


###################################################
### code chunk number 39: Ieff
###################################################
nd <- data.frame(tfI = seq(0, 15, , 151), lex.Cst = "Ins", sex = "M")
nr <- data.frame(tfI =     2            , lex.Cst = "Ins", sex = "M")
# We need to use xvars="age" in ci.exp because age is in the model
# but not in the prediction frames nd and nr
ppr <- ci.exp(pm, list(nd, nr), xvars = "age")
cpr <- ci.exp(cm, list(nd, nr))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(nd$tfI, cbind(ppr, cpr), plot = T,
         lty = c(1, 2), lwd = 3, log = "y",
         xlab = "Time since insulin (years)",
         ylab = "Mortality rate ratio")
abline(h = 1, lty = 3)


###################################################
### code chunk number 40: IeffR
###################################################
nd <- data.frame(tfI = seq(0, 15, , 151), lex.Cst = "Ins", sex = "M")
nr <- data.frame(tfI =     0            , lex.Cst = "DM" , sex = "M")
ppr <- ci.exp(pm, list(nd, nr), xvars = "age")
cpr <- ci.exp(cm, list(nd, nr))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(nd$tfI, cbind(ppr, cpr), lwd = 3,
         xlab = "Time since insulin (years)",
         ylab = "Rate ratio relative to non-Insulin",
         lty = c(1, 2), log = "y", plot = TRUE)
abline(h = 1, lty = 3)


###################################################
### code chunk number 41: 01flup.rnw:1085-1090
###################################################
ii <- glmLexis(tsNA20(dmCs),
                formula = ~ Ns(age      , knots = a.kn)
                          + Ns(      tfI, knots = i.kn)
                          + Ns(age - tfI, knots = a.kn)
                          + lex.Cst + sex)


###################################################
### code chunk number 42: 01flup.rnw:1099-1105
###################################################
im <- glmLexis(tsNA20(dmCs),
                formula = ~ Ns(age      , knots = a.kn)
                          + Ns(      tfI, knots = i.kn)
                + lex.Cst : Ns(age - tfI, knots = a.kn)
                          + lex.Cst + sex)
ci.exp(im)


###################################################
### code chunk number 43: 01flup.rnw:1114-1115
###################################################
anova(ii, im, test = 'Chisq')


###################################################
### code chunk number 44: dur-int
###################################################
pii <- ci.pred(im, ndI)
pia <- ci.pred(im, ndA)
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0) / 1.6, las = 1, bty = "n")
matshade(ndI$age, pii * 1000, plot = T, log = "y",
         xlab = "Age", ylab = "Mortality per 1000 PY",
         lty = 1, lwd = 2, col = c("blue", "forestgreen", "red"), alpha = 0.1)
matshade(ndA$age, pia * 1000)


###################################################
### code chunk number 45: dur-int-RR
###################################################
ndR <- transform(ndI, tfI = 0, lex.Cst = "DM")
cbind(head(ndI), head(ndR))
Rii <- ci.exp(im , list(ndI, ndR))
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, Rii, plot = T, log = "y",
         xlab = "Age (years)", ylab = "Rate ratio vs, non-Insulin",
         lty = 1, lwd = 2, col = c("blue", "forestgreen", "red"), alpha = 0.1)
abline(h = 1)


###################################################
### code chunk number 46: 01flup.rnw:1160-1173
###################################################
dmd <- glmLexis(dmCs,
                 from = "DM", to = "Dead",
                 formula = ~ Ns(age, knots = a.kn)
                           + sex)
ind <- glmLexis(dmCs,
                 from = "Ins", to = "Dead",
                 formula = ~ Ns(age      , knots = a.kn)
                           + Ns(      tfI, knots = i.kn)
                           + Ns(age - tfI, knots = a.kn)
                           + sex)
ini <- ci.pred(ind, ndI)
dmi <- ci.pred(dmd, ndI)
dma <- ci.pred(dmd, ndA)


###################################################
### code chunk number 47: sep-mort
###################################################
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, ini * 100, plot = TRUE, log = "y",
         xlab = "Age (years)", ylab = "Mortality rates per 100 PY",
         lwd = 2, col = "red")
matshade(ndA$age, dma*100,
         lwd = 2, col = "black")


###################################################
### code chunk number 48: sep-HR
###################################################
par(mar = c(3, 3, 1, 1), mgp = c(3, 1, 0)/1.6, las = 1, bty = "n")
matshade(ndI$age, ci.ratio(ini, dmi), plot = TRUE, log = "y",
         xlab = "Age (years)", ylab = "RR insulin vs. no insulin",
         lwd = 2, col = "red")
abline(h = 1)


###################################################
### code chunk number 49: 01flup.rnw:1222-1229
###################################################
dmCs <- cutLexis(data = dmS2,
                  cut = dmS2$doins,
            timescale = "per",
            new.state = "Ins",
            new.scale = "tfI",
         split.states = TRUE)
summary(dmCs)


###################################################
### code chunk number 50: box4
###################################################
boxes(dmCs, boxpos = list(x = c(15, 15, 85, 85),
                          y = c(85, 15, 85, 15)),
      scale.R = 1000, show.BE = TRUE)
legendbox(70, 50)


###################################################
### code chunk number 51: 01flup.rnw:1261-1267
###################################################
dmM <- mcutLexis(dmL,
           timescale = "per",
                  wh = c("doins", "dooad"),
          new.states = c("Ins", "OAD"),
          new.scales = c("tfI", "tfO"),
        ties.resolve = TRUE)


###################################################
### code chunk number 52: 01flup.rnw:1271-1274
###################################################
levels(dmM)
dmM <- Relevel(dmM,  c("DM", "OAD", "Ins", "OAD-Ins", "Ins-OAD", "Dead"))
summary(dmM, t = T)


###################################################
### code chunk number 53: 01flup.rnw:1278-1281
###################################################
wh <- c(subset(dmM, lex.Cst == "Ins-OAD")$lex.id[1:2],
        subset(dmM, lex.Cst == "OAD-Ins")$lex.id[1:2])
print(subset(dmM, lex.id %in% wh), nd = 2)


###################################################
### code chunk number 54: mbox
###################################################
boxes(dmM, boxpos = list(x = c(15, 40, 40, 85, 85, 80),
                         y = c(50, 90, 10, 90, 10, 50)),
           scale.R = 1000, show.BE = TRUE)
legendbox(6, 95)


###################################################
### code chunk number 55: mboxr
###################################################
summary(dmMr <- Relevel(dmM, list(1, 2, 3, 'OAD+Ins' = 4:5, 6)))
boxes(dmMr, boxpos = list(x = c(15, 15, 85, 85, 50),
                          y = c(85, 15, 85, 15, 50)),
            scale.R = 1000, show.BE = TRUE)


###################################################
### code chunk number 56: 01flup.rnw:1333-1342
###################################################
dmMs <- splitMulti(dmMr, age = 0:100)
summary(dmMs)
levels(dmMs)
rateIns <- gamLexis(dmMr, ~ s(age) + lex.Cst, from = 1:2   , to = 3:4)
rateOAD <- gamLexis(dmMr, ~ s(age) + lex.Cst, from = c(1,3), to = c(2, 4))
rateDth <- gamLexis(dmMr, ~ s(age) + lex.Cst)
ci.exp(rateIns, subset = "lex")
ci.exp(rateOAD, subset = "lex")
ci.exp(rateDth, subset = "lex")


###################################################
### code chunk number 57: 01flup.rnw:1511-1515
###################################################
ende <- Sys.time()
cat("  Start time:", format(anfang, "%F, %T"),
  "\n    End time:", format(  ende, "%F, %T"),
  "\nElapsed time:", round(difftime(ende, anfang, units = "mins"), 2), "minutes\n")


