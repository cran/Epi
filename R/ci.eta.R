ci.eta <- function(form, cf, vcv, # model formula, coef and vcov from fitted model
                         newdata, # prediction frame (or list of frames)
                      name.check = TRUE, # report if model matrix names as those of coef
                   alpha = 0.05, df = Inf, raw = FALSE)
{
# compute c.i. from formula and coefficients for a newdata

# allow rows of missing th the result
org.op <- options(na.action = "na.pass")
on.exit(options(org.op))

# only one-sided formula to model.matrix
if (length(form) == 3) form <- form[-2]

# matrix to multiply to the parameter vector and vcov
if (is.list(newdata))
   {
   if (!(length(newdata) %in% c(2,4)))
      stop("newdata as a list must have length 2 or 4)")
   mmx <- model.matrix(form, data = newdata[[1]])
   mmr <- model.matrix(form, data = newdata[[2]])
   mm <- mmx - mmr
   if (length(newdata) == 4)
      {
      mmx <- model.matrix(form, data = newdata[[3]])
      mmr <- model.matrix(form, data = newdata[[4]])
      mm <- mm - (mmx - mmr)
      }
   } else mm <- model.matrix(form, data = newdata)

# check sanity of formula and coeffiets
if (ncol(mm) != length(cf))
   stop("mismatch of formula and no. coefficients:\n",
        "ncol(model matrix) = ", ncol(mm),
        " but length(cf) = ", length(cf), "\n")
if (any(names(cf) != colnames(vcv)) |
    any(names(cf) != rownames(vcv)))
   stop("names of cf do not match row/col names of vcv\n")
if ((any(colnames(mm) != names(cf))) & name.check)
   {
   cat("NOTE: colnames(mm) and names(cf) do not match:\n")
   print(cbind(model = colnames(mm),
                  cf = names(cf)))
   }

# singular models produce NAs so remove from cf and vcov
out <- which(is.na(cf))
cf  <-  cf[-out]
vcv <- vcv[-out, -out]
 mm <-      mm[, -out]

# predicted eta and its vcov or ci
ct <- mm %*% cf
vc <- mm %*% vcv %*% t(mm)
se <- sqrt(diag(vc))

# confidence intervals
ci <- cbind(ct, se) %*% ci.mat(alpha = alpha, df = df)

# return results
if (raw) list(eta = ct,
              var = vc) else ci
}
