# The addCov method
addCov <- function (Lx, ...) UseMethod("addCov")

addCov.default <-
addCov.Lexis <-
function(Lx,
       clin,
  timescale = 1,
      exnam,
        tfc = "tfc")
{
# Function to add clinically measured covariates to a Lexis object

# The point is to cut the Lexis object at the examination dates
# and subsequently add the clinical records

# to avoid notes in check
org.Cst <- NULL
org.Xst <- NULL
lex.Cst <- NULL
lex.Xst <- NULL

# ...but first the usual cheking of paraphernalia

if(!inherits(Lx  ,"Lexis")) stop("Lx must be a Lexis object.\n")
if( inherits(clin,"Lexis")) stop("clin cannot be a Lexis object.\n")
    
# Is the timescale argument a timescale in Lx and is it a variable in clin? 
ts <- if(is.numeric(timescale)) timeScales(Lx)[timescale] else timescale
if(!(ts %in% timeScales(Lx)))
    stop("timescale argument (", ts,
         ") must be among the timescales in the Lexis object ",
         deparse(substitute(Lx)),":", timeScales(Lx), ".\n" )

clin.nam <- deparse(substitute(clin))
if( !( ts %in% names(clin) & "lex.id" %in% names(clin) ) )
    stop("'lex.id' and timescale '", ts,
         "' must be variables in the clin object ",
         clin.nam, "\n" )

# order clin to get the possible construction of examination names ok
clin <- clin[order(clin$lex.id, clin[,ts]),]

# check that examination dates are unique within persons
if(any(dd <- duplicated(clin[,c("lex.id",ts)])))
  {
  warning("Examination dates must be unique within persons\n",
          sum(dd), " records with duplicate dates from clin object ",
          clin.nam, " excluded.")
  clin <- clin[!dd,]
  }
    
# the variable holding the name of the examination
if(missing(exnam)) exnam <- "exnam"
# and if it is not there, construct it
if(!(exnam %in% names(clin)))
   clin[,exnam] <- paste("ex",
                          ave(clin$lex.id,
                              clin$lex.id,
                              FUN = function(x) cumsum(x/x)),
                          sep="" )
# exnam cannot have values that are also states
if( length(common <- intersect(levels(Lx),
                               unique(clin[,exnam]))) )
  stop("Levels of Lx and examination names in clin must be disjoint",
       "\nbut", paste(common, collapse=", "), "are in both")

#...done checking 
    
# variables to merge by
mvar <- c("lex.id", ts)
    
# clinical variables to be merged in
# --- note we take examination date and name as a cinical variable too 
cvar <- setdiff(names(clin), mvar)

# A data frame of cutting times of the examinations
cfr <- data.frame(lex.id = clin$lex.id,
                     cut = clin[,ts],
               new.state = clin[,exnam])

# a copy of Lx with a saved copy of the state variables in org.
Lc <- transform(Lx, org.Cst = lex.Cst,
                    org.Xst = lex.Xst)
    
# Now cut Lc at each new examination date, state variables will be
# changed to examination names
Lc <- rcutLexis(Lc,
               cut = cfr,
         timescale = ts)

# Lc now has the exnam in the variable lex.Cst, so we can merge the
# clinical data to this if we rename to exnam
# the lex.Cst contains the examination name, except where
# the original levels are left
Lc[,exnam] <- as.factor(ifelse(Lc$lex.Cst %in% levels(Lx),
                               NA,
                               as.character(Lc$lex.Cst)))
mvar <- c("lex.id", exnam)
    
# timescale is present in both Lc and clin,
# so rename in clin, it will be the date of clin 
names(clin)[grep(ts, names(clin))] <- tfc

# merge with clinical measurements keeping the attributes
att.Lc <- attributes(Lc)
Lc <- left_join(Lc, clin, by = mvar)
att.Lc$names <- attributes(Lc)$names
# attributes(Lc) <- att.Lc
    
# compute time since last examination
Lc[,tfc] <- Lc[,ts] - Lc[,tfc]

# we have problems if a clical date falls in an interval
# ending with a transition, they need to be fixed:
wh <- with(Lc, which(c(org.Xst[-nrow(Lc)] != org.Cst[-1] &
                       lex.id [-nrow(Lc)] == lex.id [-1])))
Lc$org.Xst[wh] <- Lc$org.Cst[wh+1]

# move the original states back
Lc <- select(Lc, -lex.Cst,
                 -lex.Xst) %>%
      rename(lex.Cst = org.Cst,
             lex.Xst = org.Xst)
    
# Add tfc as a time.scale, time.since and breaks:
attr(Lc, "time.scales") <- c(attr(Lx, "time.scales"), tfc) 
attr(Lc, "time.since" ) <- c(attr(Lx, "time.since" ), "" )
brt <- list(x = NULL)
names(brt) <- tfc
attr(Lc, "breaks") <- c(attr(Lx, "breaks"), brt) 
attr(Lc, "class") <- c("Lexis","data.frame")
    
# Done! - well order first
sortLexis(Lc)
}
