unLexis <-
function(Lx)
{
if (!inherits(Lx, "Lexis")) stop("Not a Lexis object")
attr(Lx, "time.scales") <- NULL
attr(Lx, "time.since") <- NULL
attr(Lx, "breaks") <- NULL
attr(Lx, "class") <- setdiff(attr(Lx, "class"), "Lexis")
Lx
}

