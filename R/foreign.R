# The msdata method
msdata <- function (obj, ...) UseMethod("msdata")

msdata.Lexis <-
function( obj,
   time.scale = timeScales(obj)[1], ... )
{
tr.mat <- tmat(obj)
# Essentially a msdata object is a stacked Lexis object with different variable names
tmp <- stack.Lexis( factorize.Lexis(obj) )
lv  <- c( match(timeScales(obj),names(tmp)),
           grep("lex\\."       ,names(tmp)) )
# The transitions that we refer to are extracted from lex.Tr:
ss <- strsplit( as.character(tmp$lex.Tr), "->" )
# The resulting dataframe is created by renaming columns in the stacked Lexis object
data.frame( id = tmp$lex.id,
          from = sapply( ss, FUN=function(x) x[1] ),
            to = sapply( ss, FUN=function(x) x[2] ),
         trans = as.integer( tmp$lex.Tr ),
        Tstart = tmp[,time.scale],
         Tstop = tmp[,time.scale] + tmp$lex.dur,
          time = tmp$lex.dur,
        status = as.integer( tmp$lex.Fail ),
                 tmp[,-lv] )
}

# The etm method
etm <- function (data, ...) UseMethod("etm")
 
etm.Lexis <-
function( data,
    time.scale = timeScales(data)[1],
     cens.name = "cens",
             s = 0,
             t = "last",
    covariance = TRUE,
      delta.na = TRUE,
           ...
          )
{
dfr <- data.frame( id = data$lex.id,
                 from = as.character(data$lex.Cst),
                   to = as.character(data$lex.Xst),
                entry = data[,time.scale],
                 exit = data[,time.scale] + data$lex.dur,
     stringsAsFactors = FALSE )
dfr$to <- with( dfr, ifelse( from==to, cens.name, to ) )
etm::etm( data = dfr,
   state.names = levels( data$lex.Cst ),
           tra = tmat(data,mode="logical"),
     cens.name = cens.name,
             s = s,
             t = t,
    covariance = covariance,
      delta.na = delta.na )
}
