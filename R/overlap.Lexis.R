# The overlap method
overlap <- function (Lx, ...) UseMethod("overlap")

overlap.default <-
overlap.Lexis <-
function(Lx, ...)
{
# must be a Lexis object
if (!inherits(Lx, "Lexis")) stop("Lx must be a Lexis object")
# the Lexis object should be sorted
Lx <- sortLexis(Lx)
# a timescale without missing values
whs <- names(timeSince(Lx))[timeSince(Lx) == ""][1]
tsc <- Lx[, timeScales(Lx)[1]]
# the next value on the timescale
tnx <- c(tsc[-1], NA)
# is this record overlapping the next and is it the same person?
olp <- ((tsc + Lx$lex.dur) > tnx) &
       duplicated(Lx$lex.id, fromLast = TRUE)
# the ids of persons with any overlapping records
if (any(olp))
   {
   who <- unique(Lx$lex.id[olp])
   cat("Overlapping records found in",
       length(who), "persons - returned invisibly.\n")
   return(invisible(who))
   }
else
   {
   cat("No overlapping records within persons.\n")
   return(NULL)
   }
}
