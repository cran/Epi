# Survfit requires a special value for censorings not just the same as
# is the convention in Lexis objects, and moreover this must be the
# first level of factor that goes into the event argument of Surv:

mkcens.Lexis <-
function(Lx, cens = "cens")
{
Rx <- Epi::sortLexis(Lx)
last <- rev(!duplicated(rev(Rx$lex.id)))
Rx$lex.Xst <- ifelse(last & Rx$lex.Cst == Rx$lex.Xst,
                     cens,
                     as.character(Rx$lex.Xst))
Relevel(factorize(Rx), cens)
}

# The AaJ method
AaJ <- function(Lx, ...) UseMethod("AaJ")

# The actual function
AaJ.default <-
AaJ.Lexis <-
function(Lx,
         formula   = ~ 1,
         timeScale = 1)
{
lex.id  <- NULL
lex.Cst <- NULL
if (!inherits(Lx, "Lexis"))
    stop("1st argument must be a Lexis object")
if( length(formula) != 2 )
    stop("'formula' must be a one-sided formula")

rhs <- as.character(formula[length(formula)]) # right hand side of formula
Lx <- mkcens.Lexis(Lx)
    
ts <- check.time.scale(Lx, timeScale)
Lx$zeit <- Lx[,ts]
cat("NOTE: Timescale is ", ts, 
#   "; initial level assumed to be ", levels(Lx)[2], 
    "\n", sep = "")

form <- Surv(Lx$zeit,
             Lx$zeit + Lx$lex.dur,
             Lx$lex.Xst) ~ 1
form[3] <- formula[2]
survfit(form, 
        id = lex.id,
    istate = lex.Cst,
      data = Lx)
}
