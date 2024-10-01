# convenience function to compute midpoints between points in a vector
midpoint <- function(x) x[-1] - diff(x) / 2

#----------------------------------------------------------------------
# Simulate rates from a list of models, 0th simulation is the rates
# based on the model estimates
simrates <-
function(mods, # named list of models - one for each cause
           nd, # prediction data frame - common for all causes
         tnam,
           nB = 1000, # number of simulated samples
          int = mean(diff(nd[,tnam]))) # time is first in column of nd
#        seed) not implemented yet
{
# from a list of models (mods) for nC cause-specific rates and a
# common prediction frame, nd, we generate nB samples of predicted
# rates. It is assumed that the nd refer to endpoints of intervals
res <- NArray(list(time = nd[,tnam],
                    mod = names(mods),
                    sim = 1:nB))
nT <- nrow(nd) - 1 # no of interval midpoints
nC <- length(mods) # no. competing risks
for(i in 1:nC) res[,i,] <- ci.lin(mods[[i]], nd, sample = nB)
res <- exp(res)    # esimated rates at the interval endpoints
attr( res, "int" ) <- int
res
}

#----------------------------------------------------------------------
# Compute nC+1 cumulative risks from nC rates
mkprobs <-
function(rtab, int)
{
# rtab is assumed to be a matrix with rates in nC (named) columns
# at nT equidistantly placed times at 0, int, 2*int, ...
# so of dimension nT x nC
nT <- nrow(rtab)
nC <- ncol(rtab)

# intensities at interval midpoints so (nT-1) x nC
mtab <- (rtab[-1,,drop=FALSE] + rtab[-dim(rtab)[1],,drop=FALSE]) / 2

# integrated intensities at interval boundaries, so nT x nC
Rtab <- rbind(0, apply(mtab, 2, cumsum)) * int

# state probabilities at interval endpoints so nT x (1+nC)
Pr <- Rtab[,c(1,1:nC)] * NA
colnames(Pr)[1] <- "Surv" # Survival prob in first column

# survival function (also the first column in Pr)
Pr[,"Surv"] <- Sv <- exp(-apply(Rtab, 1, sum))

# convenience function to compute midpoints between points in a vector
midpoint <- function(x) x[-1] - diff(x) / 2
# cumulative risks - note we are using the midvalue of the Sv, length nT
for( i in 1:nC ) Pr[,i+1] <- c(0, cumsum(mtab[,i] * midpoint(Sv)) * int)
Pr
}

#----------------------------------------------------------------------
# return first, mean and median and ci from a simulated sample
mnqt <-
function(x, alpha = 0.05)
{
# extracts the relevant quantiles across the
# simulated samples
c(quantile(x,
           probs = c(0.5, alpha / 2, 1 - alpha / 2),
           na.rm = TRUE))
}

#----------------------------------------------------------------------
# Here is the function we need to compute cumulative risks
ci.Crisk <-
function(mods, # list of models
           nd, # prediction data frame for rates at points
               # representing interval endpoints
         tnam = names(nd)[1],
           nB = 1000,   # no of parametic bootstrap samples
         perm = length(mods):0 + 1, # order of cumulation of states
        alpha = 0.05,   # 1 minus confidence level
      sim.res = 'none')
#        seed = ) not implemented yet
{
int <- mean(diff(nd[,tnam])) # interval length in nd
if (length(table(round(diff(nd[,tnam]), 6))) > 1)
   stop("Time column column of nd must be equidistant times\n")
cat("NOTE: Times are assumed to be in the column", tnam,
    "at equal distances of", int, "\n")

# First generate the simulated probabilities T x (C+1) x nB; simrt
# contains the estimates as the first entry in the simulation dimension
simrt <- simrates(mods = mods,
                    nd = nd,
                  tnam = tnam,
                    nB = nB,
                   int = int)
#                 seed = seed) not implemented yet

# if we want the simulation results for further processing
attr(simrt, "int") <- int # the interval length needed with simrt
attr(simrt, "time") <- nd[,tnam]
if(substr(tolower(sim.res), 1, 1) == 'r') return(invisible(simrt))

# Array for the state probabilities
probs <- NArray(list(time = nd[,tnam],
                    cause = c("Surv", dimnames(simrt)[["mod"]]),
                      sim = 1:nB))
names(dimnames(probs))[1] <- tnam
# Compute the cumulative risks
  # a trick to avoid dropping dimension 2 if it has length 1
  x2 <- c(1, 1:dim(simrt)[2])
for(i in 1:nB) probs[,,i] <- mkprobs(simrt[,x2,i][,-1,drop=FALSE], int = int)
attr(probs, "int") <- int # the interval length needed with probs

if(substr(tolower(sim.res), 1, 1) == 'c') return(invisible(probs))

# otherwise we compute confidence intervals for selected statistics
  # 1. state (cumulative) probabilities
  Crisk <- sim2ci.Crisk(probs, alpha = alpha)
  # 2. stacked (cumulative) probabilities
  Srisk <- sim2ci.Srisk(probs,  perm = perm,
                               alpha = alpha)
  # 3. sojourn times in states
  Stime <- sim2ci.Stime(probs, alpha = alpha)

# finally return in a list of useful quantities
rlist <- list(Crisk = Crisk,
              Srisk = Srisk,
              Stime = Stime,
               time = nd[,tnam])
attr(rlist, "int") <- int # the interval length used
return(invisible(rlist))
}

#----------------------------------------------------------------------
# 1. cumulative probabilities
sim2ci.Crisk <-
function(probs,
         alpha = 0.05)
aperm(apply(probs, 1:2, mnqt, alpha), c(2,3,1))

#----------------------------------------------------------------------
# 2. stacked (cumulative) probabilities
sim2ci.Srisk <-
function(probs,
          perm = 1:dim(probs)[2],
         alpha = 0.05)
{
cprobs <- aperm(apply(probs[,perm,],
                      c(1, 3),
                      cumsum),
                c(2,1,3))
cnam <- dimnames(cprobs)[["cause"]]
for(i in 2:length(cnam))
    dimnames(cprobs)[["cause"]][i] <- paste(cnam[1:i],
                                            collapse = "+")
cpr <- aperm(apply(cprobs, 1:2, mnqt, alpha), c(2,3,1))
cpr[,-dim(cpr)[2],]
}

#----------------------------------------------------------------------
# 3. sojourn times in states
sim2ci.Stime <-
function(probs,
           int = attr(probs, "int"),
         alpha = 0.05)
{
sojrn <- apply(probs,
                 2:3,
               function(x) cumsum(midpoint(x))) * int
res <- aperm(apply(sojrn, 1:2, mnqt, alpha), c(2,3,1))
res <- res[c(1, 1:dim(res)[1]), , ]
res[1,,] <- 0
      dimnames(res) [1] <-       dimnames(probs) [1]
names(dimnames(res))[1] <- names(dimnames(probs))[1]
res
}
