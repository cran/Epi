# Utility: which records can be merged with the next and still have
# lex.dur < lim. Assumes that variable keep.rec is in Lx
close2next <-
function(Lx, lim)
{
if (length(lim) == 1) lim <- c(lim, 3 * lim)
setdiff(
with(Lx, which(lex.Cst == lex.Xst           &
                lex.id == c(lex.id[-1], NA) &
                lex.dur < lim[1]            &
               (lex.dur + c(lex.dur[-1], NA))
                        < lim[2])),
         which(Lx$keep.rec) - 1)
}

# Utility: merge records in wh with the records in wh+1
mergeWnext <-
function(Lx, wh)
{
if (any(near(diff(wh), 1))) stop("consecutive 'wh's not allowed\n")
Lx$lex.Xst[wh] <-                  Lx$lex.Xst[wh + 1]
Lx$lex.dur[wh] <- Lx$lex.dur[wh] + Lx$lex.dur[wh + 1]
Lx[-(wh + 1),] 
}

# Utility: remove those indices that should not be used (i.e. merged
# with the next one) by removing every other among consecutive indices 
thinWh <-
function(wh)
{
dwh <- diff(c(-1, wh)) # difference between indices; any negative
                       # number will work
is1 <- near(dwh, 1)             # those with a usable successor
st1 <- near(diff(c(0, is1)), 1) # start of runs of 1
dfr <- data.frame(is1 = as.numeric(is1)) 
sq1 <- (group_by(dfr, cumsum(st1)) %>%
        mutate(z = cumsum(is1) * is1))$z
# sq1 now has runs of 1s numbered 0,1,2,3, and stand alone numbers are
# 0, so we just remove the components of wh with odd values of sq1:
wh[(sq1 %% 2) != 1]
}

# Here is the function it's all about.
coarse.Lexis <-
function(Lx, # Lexis object
        lim, # limits for proximity
       keep = FALSE) # indicator of records not combine with any previous
{
if (length(keep) != 1 &
    length(keep) != nrow(Lx)) stop("keep must have length 1 or nrow(Lx)") 
Lx$keep.rec <- keep

Lx <- sortLexis(Lx)
wh <- close2next(Lx, lim = lim)
while (length(wh > 0))
    {
wh <- thinWh(wh)        
Lx <- mergeWnext(Lx, wh)
wh <- close2next(Lx, lim = lim)
    }
Lx
}
