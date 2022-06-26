Lexis2msm <-
function(Lx,  state = "state", verbose = FALSE)
{
# the crap to pass the check
lex.dur <- lex.Xst <- lex.Cst <- NULL
    
# what are the timescales
tscal <- timeScales(Lx)
if (state %in% tscal)
   stop("The state variable cannot have the name of a time scale\n")

# find a timescale with no missing values
ts <- tscal[match(TRUE,
                  apply(Lx[,tscal],
                        2,
                        function(x) all(!is.na(x))))]

# sort so duplicated() makes sense to select last record for each person
Lx <- sortLexis(Lx)

# find the last records among consecutive records within persons
last.rec <- which(rev(!duplicated(rev(Lx$lex.id))) |
                  c( Lx[-1, ts] !=
                    (Lx[   ,ts] + Lx$lex.dur)[-nrow(Lx)],
                    TRUE))
Ll <- Lx[last.rec,]

# timescales at end of last record interval
for (i in tscal) Ll[,i] <- Ll[,i] + Ll$lex.dur

# state at start of each interval
Lx[,state] <- Lx$lex.Cst
# --- and end of last, note we carry covariates forward
Ll[,state] <- Ll$lex.Xst

# combine the data sets and remove irrelevant lex.-variables
Lmsm <- ( rbind(Ll, Lx)
      %>% sortLexis
      %>% select(-c(lex.dur, lex.Xst, lex.Cst))
        )

# re-order the columns
whLx <- c("lex.id", tscal, state)
Lmsm <- Lmsm[,c(whLx, setdiff(names(Lmsm), whLx))]

# fix the attributes
cls <- attr(Lmsm, "class")
cls[cls == "Lexis"] <- "msmLexis"
attr(Lmsm, "class") <- cls
attr(Lmsm, "breaks") <- NULL

# tell what has been going on
if (verbose) cat("Object with state occupied at time points",
                 "on the time scales:\n", tscal, "\n")
Lmsm
}
