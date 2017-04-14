\name{mcutLexis}
\alias{mcutLexis}
\title{
Cut follow-up at multiple event dates and keep track of order of events 
}
\description{
A generalization of \code{\link{cutLexis}} to the case where different
events may occur in any order.
}
\usage{
mcutLexis( L0, timescale = 1, wh,
              new.states = NULL,
        precursor.states = NULL,
              seq.states = TRUE,
              new.scales = NULL,
            ties.resolve = FALSE )
}
\arguments{
  \item{L0}{A Lexis object.}
  \item{timescale}{Which time scale do the variables in \code{L0[,wh]}
    refer to. Can be character or integer.}
  \item{wh}{Which variables contain the event dates. Character or
    integer vector}
  \item{new.states}{Names of the events forming new states. If
    \code{NULL} equal to the variable names from \code{wh}.}
  \item{precursor.states}{Which states are precursor states. See
    \code{\link{cutLexis}} for definition of precursor states.}
  \item{seq.states}{Should the sequence of events be kept track of? That
    is, should A-B be considered different from B-A. If \code{FALSE}, the
    state with combined preceding events A and B will be called A+B.}
  \item{new.scales}{Should we construct new time scales indicating the
    time since each of the event occurrences.}
  \item{ties.resolve}{Should tied event times be resolved by adding
    random noise to tied event dates. If \code{FALSE} the function will
    not accept that two events occur at the same time for a person
    (ties). If \code{TRUE} a random quantity in the range
    \code{c(-1,1)/100} will be added to all event times in all records
    with at least one tie. If numeric a random quantity in the range
    \code{c(-1,1)*ties.resolve} will be added to all event times in all
    records with at least one tie.}
}
\value{A \code{\link{Lexis}} object with extra states created by
  occurrence of a number of intermediate events.
}
\author{
Bendix Carstensen, \url{http://BendixCarstensen.com}
}
\seealso{
\code{\link{cutLexis}}, \code{\link{Lexis}}, \code{\link{splitLexis}}
}
\examples{
# A dataframe of times
set.seed(563248)
dd <- data.frame( id = 1:10,
                 doN = round(runif(10,-30, 0),1),
                 doE = round(runif(10,  0,20),1),
                 doX = round(runif(10, 50,60),1),
                 doD = round(runif(10, 50,60),1),
                 # these are the event times
                 doA = c(NA,20,NA,27,35,NA,52, 5,43,80),
                 doB = c(25,NA,37,40,NA,NA,15,23,36,61) )

# set up a Lexis object with time from entry to death/exit
Lx <- Lexis( entry = list(time=doE,
                           age=doE-doN),
              exit = list(time=pmin(doX,doD)),
       exit.status = factor(doD<doX,labels=c("OK","D")),
              data = dd )
summary( Lx )

# cut the follow-up at dates doA and doB
L2 <- mcutLexis( Lx, "time", wh=c("doA","doB"),
                 new.states = c("A","B"),
           precursor.states = "OK",
                 seq.states = TRUE,
                 new.scales = c("tfA","tfB") )
summary( L2 )
L2

# show the states
boxes( L2, boxpos=list(x=c(10,50,50,90,50,90),
                       y=c(50,90,50,90,10,10)),
           show.R=FALSE, show.BE=TRUE )


L3 <- mcutLexis( Lx, "time", wh=c("doA","doB"),
                 new.states = c("A","B"),
           precursor.states = "OK",
                 seq.states = FALSE,
                 new.scales = c("tfA","tfB") )
summary( L3 )
boxes( L3, boxpos=list(x=c(10,50,50,90,50),
                       y=c(50,90,50,50,10)),
           show.R=FALSE, show.BE=TRUE )
}
\keyword{survival}
