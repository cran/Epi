ci.eta <- function(form,    # model formula
                   cf, vcv, # coef and vcov from fitted model
                   newdata, # prediction frame
                   alpha = 0.05, df = Inf, raw = FALSE)
{
# compute c.i. from formula and coefficients for a newdata
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
if (ncol(mm) != length(coef))
   stop("mismatch of formula and coefficients:\n",
        "ncol(model matrix)=", ncol(mm),
        "length(coef)=", length(coef))
if (any(colnames(mm) != names(coef)))
   {
   cat("NOTE: colnames and coefnames do not match:\n")
   print(cbind(model = colnames(mm),
               coef = names(coef)))
   }
# predicted eta and its ci or vcov
ct <- mm %*% cf
vc <- mm %*% vcv %*% t(mm)
# confidence intervals
se <- sqrt(diag(vc))
ci <- cbind(ct, se) %*% ci.mat(alpha = alpha, df = df)
# return results
if (raw) list(eta = ct,
              var = vc) else ci
}
