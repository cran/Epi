xcutLexis <- function(Lx, cut, timescale = 1, sep = ".")
{
# the crap to pass the check
lex.Cst <- lex.Xst <-
old.Cst <- old.Xst <-
new.Cst <- new.Xst <-
 lex.id <- old.trn <- NULL
# cuts follow-up in a Lexis obejct, preserving the original states
# hence x-cut: cross-classification
Lx <- sortLexis(Lx)
Lx <- mutate(Lx, old.Cst = as.character(lex.Cst),
                 old.Xst = as.character(lex.Xst))
#
# now, cut at the new state(s), BUT
# ignoring which are absorbing states to keep new "from" states
# - absorbing states in lex.Xst are preserved in old.Xst
Lx <- cutLexis(Lx, cut,
               timescale = timescale,
               precursor.states = levels(Lx))
#
# note that old state variables have been expanded by the cutLexis, so
# if there is one or more new cuts in an initial Lx interval with a
# transition there will be multiple records with idential values of
# both old.Cst and old.Xst but old.Cst != old.Xst.
# So we need to identify such sequences of records and then change
# the value of old.Xst to that of old.Cst for these units
Lx <- mutate(Lx, new.Cst = as.character(lex.Cst),
                 new.Xst = as.character(lex.Xst),
               # different Cst and Xst and the pair identical to next
               # indicating the units to change (within person)
                 old.trn = paste0(old.Cst, old.Xst),
                 old.chg = (old.Cst != old.Xst) &
                           (lex.id  == lex.id[c(2:nrow(Lx), 1)]) &
                           (old.trn == c(old.trn[-1], "")))
Lx$old.Xst[Lx$old.chg] = Lx$old.Cst[Lx$old.chg]
#
# now for the construction of new states
Lx <- mutate(Lx, lex.Cst = ifelse(old.Cst == new.Cst,
                                  new.Cst,
                                  paste0(old.Cst, sep, new.Cst)),
                 lex.Xst = ifelse(old.Xst == new.Xst,
                                  new.Xst,
                                  paste0(old.Xst, sep, new.Xst)))
#
# claen up and beautify
todrop <- match(c("old.Cst", "new.Cst",
                  "old.Xst", "new.Xst",
                  "old.trn", "old.chg"), names(Lx))
factorize(Lx[, -todrop])
}
