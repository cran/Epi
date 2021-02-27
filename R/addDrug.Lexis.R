# The addDrug method
addDrug <- function (Lx, ...) UseMethod("addDrug")

#----------------------------------------------------------------------
# ins0 - a utility to dream up 0 purchases

ins0 <-
function(pdat, # data set of purchases with names tnam, amt, lex.id
          amt = "amt", # name of amount variable in pdat
          apt = "apt", # name of amount per time variable in pdat
                       # to be used by metod "dos"
       method = "ext", # extrapolation from times and amounts
                       # or "dos" based on dosage, apt and amt
                       # or "fix" using a fixed interval maxt
         maxt = NULL,  # maximal time covered by a single prescription,
        grace = 0,     # grace period after all used
         tnam = setdiff(names(pdat),c("lex.id", amt))[1],
      verbose = TRUE)
{
# pdat must include lex.id as variable
# inserts a 0 purchase record at end of exposure from each purchase
# 3 methods implemented:
# "ext" each prescription lasts as if consumed as the previous.
# "dos" each presription lasts amt/apt+grace.
# "fix" exch prescription lasts the same time.

# bind variables
lex.id <- dop <- extime <- NULL
    
if(length(intersect(c("lex.id", amt), names(pdat))) < 2)
  stop("lex.id and", amt, "must be columns in the data frame\n")
  
if(missing(tnam)) cat('NOTE: timescale "', tnam, '" assumed\n', sep = '')
    
if(method == "ext")
  if(verbose)
  cat("NOTE: end of exposure based on differences in purchase times (", 
      tnam,")\n and amount purchased (", amt, ").\n", sep = "")
    
if(method == "dos")
  {
  if(is.numeric(apt)) 
    {
    pdat$apt <-  apt
         apt <- "apt"
    }  
  if(!(apt %in% names(pdat))) 
    stop('method="dos" reqires amount per time, "', apt, '" in data.\n')
  if(verbose)
  cat("NOTE: end of exposure based on purchase and dosage (", apt, ").\n",
      sep = "")
  }

if(method == "fix")
  {
  if(!is.numeric(maxt)) 
    stop('method="fix" reqires a fixed coverage time in argument maxt.\n')
  if(verbose)
  cat("NOTE: end of exposure based on fixed coverage time of", maxt, ".\n")
  }

# order by purchase data within each person
pdat <- pdat[order(pdat$lex.id,pdat[,tnam]),]
    
# put the variable names in pdat for simplicity
                         pdat$dop <- pdat[,tnam]
                         pdat$amt <- pdat[, amt]
if(apt %in% names(pdat)) pdat$apt <- pdat[, apt]
    
# compute the time of expiry of current purchase by method chosen
if(method == 'fix')
  pdat$extime <- pdat$dop + maxt

if(method == 'dos')
  pdat$extime <- pdat$dop + pdat$amt / pdat$apt + grace

if(method == 'ext')
  pdat <- group_by(pdat, lex.id) %>%
          mutate(extime = dop +
                          c(NA, diff(dop) * amt[-1] /
                                            amt[-length(amt)]) 
                          + grace) %>%
          ungroup()
  
# generate 0 purchases in pzero where coverage expired before new purchase
group_by(pdat, lex.id)          %>%
filter(extime < c(dop[-1],Inf)) %>%
mutate(dop = extime, amt = 0)   %>%
select(lex.id, dop, amt)        -> pzero
    
# pzero has potentially become a grouped tibble so make it a data frame
pzero <- as.data.frame(pzero)

# append to original data frame
pdat <- rbind(pdat[,names(pzero)], pzero)

# reinstate the time-scale name
names(pdat)[grep("dop",names(pdat))] <- tnam
    
# sort by lex.id and time and coerce to data frame
oo <- order(pdat$lex.id, pdat[,tnam, drop = TRUE])
pdat <- as.data.frame(pdat[oo,])
rownames(pdat) <- NULL
pdat
}

#----------------------------------------------------------------------
# addDrug.Lexis

addDrug.default <-
addDrug.Lexis <-
function(Lx, # Lexis object, should be timesplit 
       pdat, # (named) list of data frames of drug purchaces
        amt = "amt", # name of amount variable in pdat
        apt = "apt", # name of amount per time variable in pdat
                     # to be used if method = "dos"
     method = "ext", # extrapolation from times and amounts
                     # "dos" is based on dosage, apt and amt
                     # "fix" is using a fixed interval, maxt
       maxt = NULL,  # vector of times covered by a single prescription,
      grace = 0,     # vector of grace periods after final data
       tnam = setdiff(names(pdat[[1]]),c("lex.id", amt))[1],
     prefix = NULL,  # character vector
     suffix = NULL,  # character vector
     sepfix = "."    # separator for pre- and suf-fixes
        )
{
# utility functions
  na0 <- function(x) ifelse(is.na(x), 0, x)
csum0 <- function(x) c(0, cumsum(na0(x)[-length(x)]))

# binding variables
qwzrx <- exnam <- tfc <- lex.id <- pur <-
lex.dur <- xtime <- expos <- dospt <- NULL
    
# don't bother about the warnings
oldopts <- options(warn = -1)    
on.exit(options(oldopts))
    
# Save attributes to return with result
lex.attr <- attributes(Lx)

# time scale
if (missing(tnam)) 
   cat("NOTE: timescale taken as '", tnam, "'\n", sep = "")
    
# construct names for the 0-expanded purchase data frames if not given
if(is.null(names(pdat))) names(pdat) <- paste0('P', 1:length(pdat))

# prefix or suffix for variables ex, tf, ct, cd
if(!(suff <- !is.null(prefix))) prefix <- names(pdat)    
if(!(pref <- !is.null(suffix))) suffix <- names(pdat)    
if(suff && pref) cat("NOTE: you are asking for *both* pre- and suffix\n")
# if none specified use prefix
pref <- pref | !suff

# the renaming vector to be used later 
onam <- c("expos","tfex","ctime","cdos")
shortnames <- substr(onam, 1, 2)
    
# number of purchase files 
np <- length(pdat)

# expand maxt and grace by recycling
maxt  <- rep(maxt , np)[1:np]
grace <- rep(grace, np)[1:np]
    
# Structure to hold the 0-expanded data sets
pdat0 <- list()
pall <- NULL
    
# Define the 0 amount purchases at exposure end and collect times in pall
for(i in 1:np)
   { 
   pdat0[[i]] <- ins0(pdat[[i]],        
                            amt = amt,  
                            apt = apt,  
                         method = method,
                           maxt = maxt[i], 
                          grace = grace[i],
                           tnam = tnam,
                        verbose = i == 1)
   pall <- rbind(pall, pdat0[[i]][,c("lex.id", tnam)])
   }

# order pall by id and time, remove duplicates and add a bogus variable
oo <- order(pall$lex.id, pall[,tnam])
pall <- pall[oo,]
pall <- pall[!duplicated(pall),]
pall$qwzrx <- 0   

# add the total purchase data in order to expand to all cut dates when  
# appending data from each drug
Gx <- addCov.Lexis(Lx, pall, timescale = tnam)
Gx <- select(Gx, -qwzrx, -exnam, -tfc)
    
# remove tfc as time scale and keep Lexis attributes to resinstate later
wh.tfc <- match("tfc", attr(Gx, "time.scales"))
 attr(Gx,"time.scales") <- Gsc <- attr(Gx, "time.scales")[-wh.tfc]
 attr(Gx,"time.since")  <- Gsi <- attr(Gx, "time.since" )[-wh.tfc]      
class(Gx)               <- Gcl <- c("Lexis", "data.frame")
    
# put the expanded drug purchases in the Lexis object 
allnam <- NULL
for(i in 1:np)
   {
   # renaming the variables with pre- and/or suffix
            rnam <- shortnames
   if(suff) rnam <- paste(rnam, suffix[i], sep = sepfix)
   if(pref) rnam <- paste(prefix[i], rnam, sep = sepfix)
   names(rnam) <- onam
   allnam <- c(allnam, rnam)
       
   # add the i'th drug exposure data file
    Gx <- addCov.Lexis(Gx, pdat0[[i]], timescale = tnam, exnam = "pur")
     
   # compute the exposure variables, drop unneeded variables and rename
   Gx <- 
     ( group_by(Gx, lex.id, pur)             
   %>% mutate(xtime = sum(lex.dur) + tfc[1],
                      # total exposure time for this purchase 
              dospt = amt / xtime ) # dose per time for this interval
   %>% group_by(lex.id) 
   %>% mutate(expos = !is.na(amt) & amt>0,
               tfex = csum0(lex.dur * (cumsum(expos) > 0)),
              ctime = csum0(lex.dur * expos * (amt>0)) + na0(tfc[1]),
               cdos = csum0(lex.dur * expos * dospt)
                                       + na0(dospt[1]) * na0(tfc[1]))
   %>% select(-pur, -xtime, -dospt, -tfc, -amt)
   %>% plyr::rename(rnam)
     )

   # reinstate as Lexis object       
    attr(Gx,"time.scales") <- Gsc
    attr(Gx,"time.since")  <- Gsi
   class(Gx)               <- Gcl
   }
    
return(Gx[,c(setdiff(names(Gx), allnam), allnam)])
}
