### R code from vignette source 'addLexis.rnw'

###################################################
### code chunk number 1: addLexis.rnw:20-31
###################################################
options(width = 90,
        SweaveHooks = list(fig = function()
        par(mar = c(3,3,1,1),
            mgp = c(3,1,0) / 1.6,
            las = 1,
           lend = "butt",
            bty = "n")))
library(Epi)
library(popEpi)
library(dplyr)
library(tidyr)


###################################################
### code chunk number 2: addLexis.rnw:35-37
###################################################
anfang <- Sys.time()
cat("Start time:", format(anfang, "%F, %T"), "\n")


###################################################
### code chunk number 3: addLexis.rnw:105-117
###################################################
xcoh <- structure(list(id = c("A", "B", "C"),
                    birth = c("1952-07-14", "1954-04-01", "1987-06-10"),
                    entry = c("1965-08-04", "1972-09-08", "1991-12-23"),
                     exit = c("1997-06-27", "1995-05-23", "1998-07-24"),
                     fail = c(1, 0, 1) ),
                   .Names = c("id", "birth", "entry", "exit", "fail"),
                row.names = c("1", "2", "3"),
                    class = "data.frame" )
xcoh$dob <- cal.yr(xcoh$birth)
xcoh$doe <- cal.yr(xcoh$entry)
xcoh$dox <- cal.yr(xcoh$exit )
xcoh


###################################################
### code chunk number 4: addLexis.rnw:124-132
###################################################
Lcoh <- Lexis(entry = list(per = doe),
               exit = list(per = dox,
                           age = dox - dob),
                 id = id,
        exit.status = factor(fail, 0:1, c("Alive","Dead")),
               data = xcoh)
str(Lcoh)
(Lx <- Lcoh[,1:6])


###################################################
### code chunk number 5: addLexis.rnw:148-151
###################################################
Lx$lex.id <- as.character(Lx$lex.id)
str(Lx)
Lx


###################################################
### code chunk number 6: addLexis.rnw:159-171
###################################################
clin <- data.frame(lex.id = c("A", "A", "C", "B", "C"),
                      per = cal.yr(c("1977-3-17",
                                     "1973-7-29",
                                     "1996-3-1",
                                     "1990-7-14",
                                     "1989-1-31")),
                       bp = c(120, 140, 160, 157, 145),
                     chol = c(NA, 5, 8, 9, 6),
                     xnam = c("X2", "X1", "X1", "X2", "X0"),
         stringsAsFactors = FALSE)
str(clin)
clin


###################################################
### code chunk number 7: addLexis.rnw:185-186
###################################################
(Cx <- addCov.Lexis(Lx, clin))


###################################################
### code chunk number 8: addLexis.rnw:200-202
###################################################
(Dx <- addCov.Lexis(Lx, clin, exnam = "xnam", tfc = "tfCl"))
summary(Dx, t=T)


###################################################
### code chunk number 9: addLexis.rnw:209-224
###################################################
# split BEFORE add
Lb <- addCov.Lexis(splitLexis(Lx,
                      time.scale = "age",
                          breaks = seq(0, 80, 5)),
                   clin,
                   exnam = "xnam" )
Lb
#
# split AFTER add
La <- splitLexis(addCov.Lexis(Lx,
                            clin,
                           exnam = "xnam" ),
                 time.scale = "age",
                     breaks = seq(0, 80, 5))
La


###################################################
### code chunk number 10: addLexis.rnw:230-233
###################################################
La$tfc == Lb$tfc
La$age == Lb$age
La$per == Lb$per


###################################################
### code chunk number 11: addLexis.rnw:236-248
###################################################
## split BEFORE add
Mb <- addCov.Lexis(splitMulti(Lx, age = seq(0, 80, 5)),
                   clin,
                   exnam = "xnam" )
##
## split AFTER add
Ma <- splitMulti(addCov.Lexis(Lx,
                              clin,
                              exnam = "xnam" ),
                 age = seq(0, 80, 5))
La$tfc == Mb$tfc
Ma$tfc == Mb$tfc


###################################################
### code chunk number 12: addLexis.rnw:263-269
###################################################
cov <- c("bp", "chol")
Lx <- La
Lx <- group_by(Lx, lex.id) %>%
         fill(all_of(cov)) %>%
      ungroup()
class(Lx)


###################################################
### code chunk number 13: addLexis.rnw:274-280
###################################################
Lx <- La
Lx[,cov] <- as.data.frame(group_by(Lx, lex.id)
                          %>% fill(all_of(cov)))[,cov]
class(Lx)
La
Lx


###################################################
### code chunk number 14: addLexis.rnw:348-361
###################################################
fu <- data.frame(doe = c(2006, 2008),
                 dox = c(2015, 2018),
                 dob = c(1950, 1951),
                 xst = factor(c("A","D")))
Lx <- Lexis(entry = list(per = doe,
                         age = doe- dob),
             exit = list(per = dox),
      exit.status = xst,
             data = fu)
Lx <- subset(Lx, select = -c(doe, dob, dox, xst))
Sx <- splitLexis(Lx, "per", breaks = seq(1990, 2020, 0.6))
summary(Sx)
str(Sx)


###################################################
### code chunk number 15: addLexis.rnw:370-380
###################################################
set.seed(1952)
rf <- data.frame(per = c(2005 + runif(12, 0, 10)),
                 amt = sample(2:4, 12, replace = TRUE),
              lex.id = sample(1:2, 12, replace = TRUE)) %>%
      arrange(lex.id, per)

rg <- data.frame(per = c(2009 + runif(10, 0, 10)),
                 amt = sample(round(2:4/3,1), 10, replace = TRUE),
              lex.id = sample(1:2, 10, replace = TRUE)) %>%
      arrange(lex.id, per)


###################################################
### code chunk number 16: addLexis.rnw:392-395
###################################################
pdat <- list(F = rf, G = rg)
pdat
Lx


###################################################
### code chunk number 17: addLexis.rnw:404-417
###################################################
summary(Sx) ; names(Sx)
ex1 <- addDrug.Lexis(Sx, pdat, method = "ext") # default
summary(ex1) ; names(ex1)
print(ex1, nd = 2)
ex2 <- addDrug.Lexis(Sx, pdat, method = "ext", grace = 0.5)
summary(ex2)
print(ex2, nd = 2)
dos <- addDrug.Lexis(Sx, pdat, method = "dos", dpt = 6)
summary(dos)
print(dos, nd = 2)
fix <- addDrug.Lexis(Sx, pdat, method = "fix", maxt = 1)
summary(fix)
print(fix, nd = 2)


###################################################
### code chunk number 18: addLexis.rnw:426-435
###################################################
data(DMlate) ; str(DMlate)
Lx <- Lexis(entry = list(per = dodm,
                         age = dodm - dobth,
                         tfd = 0),
             exit = list(per = dox),
      exit.status = factor(!is.na(dodth),
                           labels = c("DM", "Dead")),
             data = DMlate[sample(1:nrow(DMlate), 1000),])
summary(Lx)


###################################################
### code chunk number 19: addLexis.rnw:439-442
###################################################
Sx <- splitLexis(Lx[,1:7], time.scale="age", breaks = 0:120)
summary(Sx)
str(Sx)


###################################################
### code chunk number 20: addLexis.rnw:450-491
###################################################
set.seed(1952)

purA <-
  ( data.frame(lex.id = rep(Lx$lex.id,
                            round(runif(nrow(Lx), 0, 20))))
%>% left_join(Lx[,c("lex.id", "dodm", "dox")])
%>% mutate(per = dodm + runif(length(dodm), -0.1, 0.99) * (dox - dodm),
           amt = sample(4:20*10, length(dodm), replace = TRUE),
           dpt = amt * round(runif(length(dodm), 3, 7)))
%>% select(-dodm, -dox)
%>% arrange(lex.id, per)
  )
addmargins(table(table(purA$lex.id)))
str(purA)

purB <-
  ( data.frame(lex.id = rep(Lx$lex.id,
                            round(pmax(runif(nrow(Lx), -10, 15), 0))))
%>% left_join(Lx[,c("lex.id", "dodm", "dox")])
%>% mutate(per = dodm + runif(length(dodm), -0.1, 0.99) * (dox - dodm),
           amt = sample(4:20*10, length(dodm), replace = TRUE),
           dpt = amt * round(runif(length(dodm), 5, 9)))
%>% select(-dodm, -dox)
%>% arrange(lex.id, per)
  ) -> purB
addmargins(table(table(purB$lex.id)))
str(purB)

purC <-
  ( data.frame(lex.id = rep(Lx$lex.id,
                            round(pmax(runif(nrow(Lx), -5, 12), 0))))
%>% left_join(Lx[,c("lex.id", "dodm", "dox")])
%>% mutate(per = dodm + runif(length(dodm), -0.1, 0.99) * (dox - dodm),
           amt = sample(4:20*10, length(dodm), replace = TRUE),
           dpt = amt * round(runif(length(dodm), 5, 7)))
%>% select(-dodm, -dox)
%>% arrange(lex.id, per)
  )
addmargins(table(table(purC$lex.id)))
str(purC)
head(purC)


###################################################
### code chunk number 21: addLexis.rnw:506-513
###################################################
Sx1 <- subset(Sx, lex.id < 100)
pur <- list(A = subset(purA, lex.id < 1000),
            B = subset(purB, lex.id < 1000),
            C = subset(purC, lex.id < 1000))
system.time(ad1 <- addDrug.Lexis(Sx1, pur, tnam = "per", grace = 1/4))
summary(Sx1)
summary(ad1)


###################################################
### code chunk number 22: addLexis.rnw:517-524
###################################################
Sx2 <- subset(Sx, lex.id < 500)
pur <- list(A = subset(purA, lex.id < 500),
            B = subset(purB, lex.id < 500),
            C = subset(purC, lex.id < 500))
system.time(ad2 <- addDrug.Lexis(Sx2, pur, tnam = "per", grace = 1/6))
summary(Sx2)
summary(ad2)


###################################################
### code chunk number 23: addLexis.rnw:531-538
###################################################
pur <- list(A = subset(purA, lex.id < 100 & runif(nrow(purA)) < 0.5),
            B = subset(purB, lex.id < 100 & runif(nrow(purB)) < 0.5),
            C = subset(purC, lex.id < 100 & runif(nrow(purC)) < 0.5))
sapply(pur, nrow)
system.time(ad3 <- addDrug.Lexis(Sx1, pur, tnam = "per", grace = 1/6))
summary(Sx1)
summary(ad3)


###################################################
### code chunk number 24: addLexis.rnw:554-560
###################################################
pur <- list(B = subset(purB, lex.id < 100),
            C = subset(purC, lex.id < 100))
sapply(pur, nrow)
system.time(ad4 <- addDrug.Lexis(Sx1, pur, tnam = "per", grace = 1/6))
summary(Sx1)
summary(ad4)


###################################################
### code chunk number 25: addLexis.rnw:569-570
###################################################
summary(ad1$lex.dur)


###################################################
### code chunk number 26: addLexis.rnw:587-590
###################################################
summary(ad1)
summary(adc <- coarse.Lexis(ad1, lim = c(1/6,1/2)))
summary(adc$lex.dur)


###################################################
### code chunk number 27: addLexis.rnw:603-613
###################################################
summary(Sx2)
system.time(ad4 <- addDrug.Lexis(Sx2,
                                 pur,
                                tnam = "per",
                               grace = 1/6))
summary(ad4)
#
ad5 <- coarse.Lexis(ad4,
                    lim = c(1/4, 1/2))
summary(ad5)


###################################################
### code chunk number 28: addLexis.rnw:618-624
###################################################
ad4$keep <- with(ad4, (B.ex & B.ct == 0) |
                      (C.ex & C.ct == 0))
ad6 <- coarse.Lexis(ad4,
                    lim = c(1/4, 1/2),
                   keep = ad4$keep)
summary(ad6)


###################################################
### code chunk number 29: addLexis.rnw:653-657
###################################################
ende <- Sys.time()
cat("  Start time:", format(anfang, "%F, %T"), "\n")
cat("    End time:", format(  ende, "%F, %T"), "\n")
cat("Elapsed time:", round(difftime(ende, anfang, units = "mins"), 2), "minutes\n")


