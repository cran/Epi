# The paths method
paths <- function (Lx, ...) UseMethod("paths")

paths.default <-
paths.Lexis <-
function(Lx, dfr = FALSE, ...)
{
# sort the Lx object
Lx <- sortLexis(Lx)
# indicator of last record for each person
last <- !duplicated(Lx$lex.id, fromLast = TRUE)
# records with transitions
tran <- with(Lx, lex.Cst != lex.Xst)
# either of the two - not as Lexis
base <- unLexis(Lx[tran | last, c("lex.id", "lex.Cst","lex.Xst")])
# convert to character
for (i in 2:3) base[,i] <- as.character(base[,i])
# split by persons
lbas <- split(base, base$lex.id)
# construct character vector of visited states, weed out duplicates
visits <- function(x)
          {
          zz <- c(x$lex.Cst[1], x$lex.Xst)
          zz[c(TRUE, zz[-1] != zz[-length(zz)])]
          }
# create list of character vectors
paths <- lapply(lbas, visits)
# paste them together
ff <- factor(sapply(paths, paste, collapse = "->"))

# if required, render as data frame else just return factor ff
if (!dfr) return(ff)
else {
     df <- data.frame(lex.id = names(ff), path = ff)
     if (is.numeric(Lx$lex.id)) df$lex.id <- as.numeric(df$lex.id)
     return(df)
     }
}
