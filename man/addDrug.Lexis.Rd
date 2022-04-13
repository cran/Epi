\name{addDrug.Lexis}
\alias{addDrug.Lexis}
\alias{coarse.Lexis}
\title{
Expand a Lexis object with information of drug exposure based on
purchase dates and -amounts
}
\description{  
A \code{\link{Lexis}} object will contain information on follow-up for a
cohort of persons through time, each record containing information of
one time interval, including the time at the beginning of each
interval. If information on drug purchase is known for the persons via
\code{lex.id} in a list of data frames, \code{addDrug.Lexis} will expand
the \code{Lexis} object by cutting at all drug purchase dates, and
compute the exposure status for any number of drugs, and add these as
variables.

In some circumstances the result is a Lexis object with a very large
number of very small follow-up intervals. The function
\code{coarse.Lexis} combines consecutive follow-up intervals using the
covariates from the first

}

\usage{
addDrug.Lexis(Lx,  # Lexis object
            pdat,  # list of data frames with drug purchase information
             amt = "amt", # name of the variable with purchased amount
             dpt = "dpt", # name of the variable with amount consumed per time
             apt = NULL,  # old name for dpt
          method = "ext", # method use to compute exposure 
            maxt = NULL,  # max duration for a purchase when using "fix"
           grace = 0,     # grace  period to be added
            tnam = setdiff(names(pdat[[1]]), c("lex.id", amt))[1],
                          # name of the time variable from Lx
          prefix = TRUE,  # should drug names prefix variable names 
          sepfix = ".",   # what should the separator be when forming prefix/suffix
         verbose = TRUE)
coarse.Lexis(Lx, lim, keep = FALSE)
}

\arguments{
  \item{Lx}{A \code{Lexis} object.
}
  \item{pdat}{Named list of data frames with drug \code{p}urchase \code{dat}a.
}
  \item{amt}{Name of the variable in the data frames in \code{pdat} with
    the purchased \code{am}oun\code{t}.
}
  \item{dpt}{Name of the variable in the data frames in \code{pdat} with
    the consumed \code{d}ose \code{p}er \code{t}ime. Must be given in
    units of units of \code{amt} per units of \code{lex.dur} in \code{Lx}.
}
  \item{apt}{Name previously used for \code{dpt}. Will disappear in next
    version.
}
  \item{method}{Character. One of \code{"ext"} (default), \code{"dos"}
    or  \code{"fix"}, for a description, see details.
}
  \item{maxt}{Numerical. Maximal duration for a purchase when using
    \code{method="fix"}, same units as \code{lex.dur}.
}
  \item{grace}{Numeric. Grace period to be added after last time of
    computed drug coverage to define end of exposure, same units as
    \code{lex.dur}. 
}
  \item{tnam}{Character. Name of the timescale used in the data frames
  in \code{pdat}. 
}
  \item{prefix}{Logical. Should the names of \code{pdat} be used as
  prefix for the 4 generated exposure variables for each drug. If false
  the names of \code{pdat} will be used as suffix.
}
  \item{sepfix}{Character, used to separate the \code{prefix} and the
    name of the generated type of variable.
}
  \item{verbose}{Logical. Should the function tell you about the choices
  you made?
}
 \item{lim}{Numeric vector of length 2. Consecutive follow-up intervals
  are combined if the first has \code{lex.dur} < \code{lim[1]}, and the
  sum of \code{lex.dur} in the two intervals is smaller than
  \code{lim[2]}. If a scalar i given, \code{c(lim,3*lim)} is used.
}
  \item{keep}{Logical of length 1 or \code{nrow(Lx)} that points to
  records that cannot be combined with preceding records.
}
}

\details{
This function internally uses \code{\link{addCov.Lexis}} to attach
exposure status for several drugs (dispensed medicine) to follow-up in a
\code{Lexis} object. Once that is done, the exposure measures are
calculated at each time.

There is one input data frame per type of drug, each with variables
\code{lex.id}, \code{amt}, a timescale variable and possibly a variable
\code{dpt}.

Three different methods for computing drug exposures from dates and
amounts of purchases are supported via the argument \code{method}.

\itemize{ 
 \item \code{"ext"}: Extrapolation: the first drug purchase is assumed
 consumed over the interval to the second purchase. Exposure
 for subsequent purchases are assumed to last as long as it would have if
 consumed at a speed corresponding to the previous purchase being
 consumed over the time span between the previous and current purchase,
 plus a period of length \code{grace}.
 
 \item \code{"dos"}: Dosage: assumes that each purchase lasts
 \code{amt}/\code{dpt} plus \code{grace}.
 
 \item \code{"fix"}: Fixed time: assumes that each purchase lasts
 \code{maxt}.
} 

So for each purchase we have defined an end of coverage (expiry
date). If next purchase is before this, we assume that the amount
purchased is consumed over the period between the two purchases,
otherwise over the period to the end of coverage. So the only difference
between the methods is the determination of the coverage for each
purchase.

Based on this, for each date in the resulting \code{\link{Lexis}} four
exposure variables are computed, see next section.
}

\value{  
A \code{\link{Lexis}} object with the same risk time, states and events
as \code{Lx}. The follow-up for each person has been cut at the purchase
times of each of the drugs, as well as at the expiry times for each drug
coverage. Further, for each drug (i.e. the data frame in the \code{pdat}
list) the name of the \code{pdat} component determines the prefix for
the 4 variables that will be added. Supposing this is \code{AA} for a
given drug, then 4 new variables will be:
\itemize{
 \item \code{AA.ex}: logical; is the person exposed in this interval
 \item \code{AA.tf}: numeric: time since first purchase, same units as
 \code{lex.dur} 
 \item \code{AA.ct}: numeric: cumulative time on the drug, same units
 as \code{lex.dur} 
 \item \code{AA.cd}: numeric: cumulative dose of the drug, same units
 as \code{amt} 
 }
So if \code{pdat} is a list of length 3 with names \code{c("a","b","c")}
the function will add variables
\code{a.ex, a.tf, a.ct, a.cd,
      b.ex, b.tf, b.ct, b.cd,
      c.ex, c.tf, c.ct, c.cd}
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link{gen.exp}},
\code{\link{addCov.Lexis}},
\code{\link{cutLexis}},
\code{\link{rcutLexis}},
\code{\link{mcutLexis}}
}
\examples{
# Follow-up of 2 persons
clear()
fu <- data.frame(doe = c(2006, 2008),
                 dox = c(2015, 2018),
                 dob = c(1950, 1951),
                 xst = factor(c("A","D")))
Lx <- Lexis(entry = list(per = doe,
                         age = doe- dob),
             exit = list(per = dox),
      exit.status = xst,
             data = fu)
Lx <- subset(Lx, select = -c(doe, dob, dox, xst))

# split FU in 1 year intervals
Sx <- splitLexis(Lx, "per", breaks = seq(1990, 2020, 1.0))

# drug purchases, one data frame for each drug 
ra <- data.frame(per = c(2007 + runif(12,0,10)),
                 amt = sample(2:4, 12, r = TRUE),
              lex.id = sample(1:2, 12, r = TRUE))
ra <- ra[order(ra$lex.id, ra$per),]

rb <- data.frame(per = c(2009 + runif(10, 0, 10)),
                 amt = sample(round(2:4/3,1), 10, r = TRUE),
              lex.id = sample(1:2, 10, r = TRUE))
rb <- rb[order(rb$lex.id, rb$per),]

# put in a named list
pdat <- list(A = ra, B = rb)
pdat

ex1 <- addDrug.Lexis(Sx, pdat, method = "ext") # default
summary(ex1)
# collapsing some of the smaller intervals with the next
summary(coarse.Lexis(ex1, c(0.2,0.5)))

ex2 <- addDrug.Lexis(Sx, pdat, method = "ext", grace = 0.2)
dos <- addDrug.Lexis(Sx, pdat, method = "dos", dpt = 6)
fix <- addDrug.Lexis(Sx, pdat, method = "fix", maxt = 1)
}
\keyword{survival}
\keyword{manip}
