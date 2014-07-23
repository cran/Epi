\name{crr.Lexis}
\alias{crr.Lexis}
\title{Fit a competing risks regression model (Fine-Gray model) using a
  Lexis object)
}

\description{
Fits a competing risks regression model using a \code{\link{Lexis}}
object assuming that every person enters at time 0 and exits at time
\code{lex.dur}. Thus is only meaningful for Lexis objects with one record
per person, (so far).
}

\usage{
crr.Lexis( obj, mod, quiet=FALSE, ...)
}

\arguments{
  \item{obj}{A Lexis object; variables in \code{mod} are taken from this.}
  \item{mod}{Formula, with the l.h.s. a character constant equal to a
    level of \code{obj$lex.Xst}, and the r.h.s. a model formula
    interpreted in \code{obj}.}
  \item{quiet}{Logical indicating whether a brief summary should be printed.}
  \item{\dots}{Further arguments passed on to \code{\link[cmprsk:crr]{crr}}.}
}

\details{
This function is a simple wrapper for \code{crr}, allowing a
formula-specification of the model (which allows specifications of
covariates on the fly), and utilizing the structure of Lexis
objects to simplify specification of the outcome. Prints a summary of
the levels used as event, competing events and censoring.

By the structure of the \code{\link{Lexis}} object it is not necessary
to indicate what the censoring code or competing events are, that is
automatically derived from the \code{Lexis} object.

Currently only one state is allowed as l.h.s. (response) in \code{mod}.
}

\value{
A \code{\link[cmprsk:crr]{crr}} object (which is a list), with two extra
elements in the list, \code{model.Lexis} - the model formula supplied,
and \code{transitions} - a table of transitions and censorings showing
which transition was analysed and which were taken as competing events.
}

\author{ Bendix Carstensen, \url{BendixCarstensen.com} }

\seealso{
  \code{\link[cmprsk:crr]{crr}}, \code{\link{Lexis}}
}

\examples{
# Thorotrats patients, different histological types of liver cancer
# Load thorotrast data, and restrict to exposed
data(thoro)
tht <- thoro[thoro$contrast==1,]
# Define exitdate as the date of livercancer
tht$dox <- pmin( tht$liverdat, tht$exitdat, na.rm=TRUE )
tht <- subset( tht, dox > injecdat )
# Convert to calendar years in dates
tht <- cal.yr( tht )

# Set up a Lexis object with three subtypes of liver cancer and death
tht.L <- Lexis( entry = list( per = injecdat,
                              tfi = 0 ),
                 exit = list( per = dox ),
          exit.status = factor( 1*hepcc+2*chola+3*hmang+
                                4*(hepcc+chola+hmang==0 & exitstat==1),
                                labels=c("No cancer","hepcc","chola","hmang","Dead") ),
                 data = tht )
summary( tht.L )

# Show the transitions
boxes( tht.L, boxpos=list(x=c(20,rep(80,3),30),
                          y=c(60,90,60,30,10) ),
              show.BE=TRUE, scale.R=1000 )

# Fit a model for the Hepatocellular Carcinoma as outcome
# - note that you can create a variable on the fly:
library( cmprsk )
hepcc <- crr.Lexis( tht.L, "hepcc" ~ volume + I(injecdat-1940) )
hepcc$model.Lexis
hepcc$transitions

# Models for the three other outcomes:
chola <- crr.Lexis( tht.L, "chola" ~ volume + I(injecdat-1940) )
hmang <- crr.Lexis( tht.L, "hmang" ~ volume + I(injecdat-1940) )
dead  <- crr.Lexis( tht.L, "Dead"  ~ volume + I(injecdat-1940) )

# Compare the effects
# NOTE: This is not necessarily a joint model for all transitions.
zz <- rbind( ci.exp(hepcc),
             ci.exp(chola),
             ci.exp(hmang),
             ci.exp(dead) )
zz <- cbind( zz[c(1,3,5,7)  ,],
             zz[c(1,3,5,7)+1,] )
rownames( zz ) <- c("hepcc","chola","hmang","dead")
colnames( zz )[c(1,4)] <- rownames( ci.exp(chola) )
round( zz, 3 )
}

\keyword{survival}
