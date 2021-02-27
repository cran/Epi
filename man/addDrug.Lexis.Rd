\name{addDrug.Lexis}
\alias{addDrug.Lexis}
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
}

\usage{
addDrug.Lexis(Lx,  # Lexis object
            pdat,  # list of data frames with drug purchase information
             amt = "amt", # name of the variable with purchased amount
             apt = "apt", # name of the variable with amount consumed per time
          method = "ext", # method use to compute exposure 
            maxt = NULL,  # max duration for a purchase when using \code{method="fix"}
           grace = 0,     # grace  period to be added
            tnam = setdiff(names(pdat[[1]]), c("lex.id", amt))[1],
                          # name of the 
          prefix = NULL,  # character vector of length \code{length(pdat)}
          suffix = NULL,  # character vector of length \code{length(pdat)}
          sepfix = ".")   # what should the separator be when forming prefix/suffix
}

\arguments{
  \item{Lx}{A \code{Lexis} object.
}
  \item{pdat}{Named list of data frames with drug \code{p}urchase \code{dat}a.
}
  \item{amt}{Name of the variable in the data frames in \code{pdat} with
    the purchased \code{am}oun\code{t}.
}
  \item{apt}{Name of the variable in the data frames in \code{pdat} with
    the consumed \code{a}mount \code{p}er \code{t}ime. Must be given in
    units of units of \code{amt} per units of \code{lex.dur} in \code{Lx}.
}
  \item{method}{Character. One of \code{"ext"} (default), \code{"amt"}
    or  \code{"fix"}, for a description, see details.
}
  \item{maxt}{Numerical. Maximal duration for a purchase when using
    \code{method="fix"}, same units as \code{lex.dur}.
}
  \item{grace}{Numeric. Grace period to be added after last time of computed drug
    coverage to define end of exposure, same units as \code{lex.dur}.
}
  \item{tnam}{Character. Name of the timescale used in the data frames in \code{pdat}.
}
  \item{prefix}{Character vector of length \code{length(pdat)}. Used to
      prefix the 4 generated exposure variables.
}
  \item{suffix}{Character vector of length \code{length(pdat)}. Used to
      suffix the 4 generated exposure variables.
}
  \item{sepfix}{Character, used to separate the \code{prefix} and the
    name of the generated type of variable.
}
}

\details{
This function internally uses \code{\link{addCov.Lexis}} to attach
exposure status for several drugs (dispensed medicine) to follow-up in a
\code{Lexis} object. Once that is done, the exposure measures are
calculated at each time.

There is one input data frame per type of drug, each with variables
\code{lex.id}, \code{amt}, a timescale variable and possibly \code{apt}.

Three different methods for computing drug exposures from dates and
amounts of purchases are supported via the argument \code{method}.

\itemize{ 
 \item \code{"ext"}: Extrapolation: the first drug purchase is assumed
 consumed over the interval to the second purchase. Exposure
 for next purchases are assumed to last as long as it would have if
 consumed at a speed corresponding to the previous purchase being
 consumed over the time span between the previous and current purchase,
 plus a period of length \code{grace}.
 
 \item \code{"dos"}: Dosage: assumes that each purchase lasts
 \code{amt}/\code{apt} plus \code{grace}.
 
 \item \code{"fix"}: Fixed time: assumes that each purchase lasts
 \code{maxt}.
} 

So for each purchase we have defined an end of coverage (expiry
date). If next purchase is before this, we assume that the amount
purchased is consumed over the period between the two purchases,
otherwise over the perido to the end of coverage. So the only difference
between the methods is the determination of the coverage for each
purchase.

Based on this, for each date in the resulting \code{\link{Lexis}} four
exposure variables are computed, see next section.
}

\value{  
A \code{\link{Lexis}} object with the same risk time, states and events
as \code{Lx}. The follow-up for each person has been cut at the purchase
times of each of the drugs, as well as at the expiry times for drug
coverage. Further for each drug (i.e. data frame in the \code{pdat}
list) 4 variables have been added:

\itemize{
 \item \code{pre.ex}: logical; is the person exposed in this interval
 \item \code{pre.tf}: numeric: time since first purchase, same units as
 \code{lex.dur} 
 \item \code{pre.ct}: numeric: cumulative time on the drug, same units
 as \code{lex.dur} 
 \item \code{pre.cd}: numeric: cumulative dose of the drug, same units
 as \code{amt} 
 }
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link{addCov.Lexis}},
\code{\link{gen.exp}},
\code{\link{cutLexis}},
\code{\link{rcutLexis}},
\code{\link{mcutLexis}}
}
\examples{
# Follow-up of 2 persons
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

rc <- data.frame(per = c(2007 + runif(12,0,10)),
                 amt = sample(2:4, 12, r = TRUE),
              lex.id = sample(1:2, 12, r = TRUE))
rc <- rc[order(rc$lex.id, rc$per),]

rd <- data.frame(per = c(2009 + runif(10, 0, 10)),
                 amt = sample(round(2:4/3,1), 10, r = TRUE),
              lex.id = sample(1:2, 10, r = TRUE))
rd <- rd[order(rd$lex.id, rd$per),]

pdat <- list(Joe = rc, Pat = rd)
pdat

addDrug.Lexis(Sx, pdat, method = "ext") # default
addDrug.Lexis(Sx, pdat, method = "ext", grace = 0.2)
addDrug.Lexis(Sx, pdat, method = "dos", apt = 6)
addDrug.Lexis(Sx, pdat, method = "fix", maxt = 1)
}
\keyword{survival}
\keyword{manip}
