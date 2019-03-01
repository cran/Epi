\name{cutLexis}
\alias{cutLexis}
\alias{countLexis}
\title{
  Cut follow-up at a specified date for each person.
}
\description{

  Follow-up intervals in a Lexis object are divided into two
  sub-intervals: one before and one after an intermediate event.  The
  intermediate event may denote a change of state, in which case the
  entry and exit status variables in the split Lexis object are
  modified.

}
\usage{
cutLexis( data, cut, timescale = 1,
                     new.state = nlevels(data$lex.Cst)+1,
                     new.scale = FALSE,
                  split.states = FALSE,
                   progressive = FALSE,
              precursor.states = NULL,
                         count = FALSE )
countLexis( data, cut, timescale = 1 )
}
\arguments{
  \item{data}{A \code{Lexis} object.}
  \item{cut}{A numeric vector with the times of the intermediate event.
    If a time is missing (\code{NA}) then the event is assumed to occur
    at time \code{Inf}. \code{cut} can also be a dataframe, see details.}
  \item{timescale}{The timescale that \code{cut} refers to. Numeric or character.}
  \item{new.state}{The state to which a transition occur at time
    \code{cut}. It may be a single value, which is then applied to all
    rows of \code{data}, or a vector with a separate value for each row}
  \item{new.scale}{Name of the timescale defined as "time since entry to
    new.state". If \code{TRUE} a name for the new scale is constructed.
    See details.}
  \item{split.states}{Should states that are not precursor states be split
    according to whether the intermediate event has occurred.}
  \item{progressive}{a logical flag that determines the changes to exit
    status. See details.}
  \item{precursor.states}{an optional vector of states to be considered
    as "less severe" than \code{new.state}. See Details below}
  \item{count}{logical indicating whether the \code{countLexis} options should
    be used. Specifying \code{count=TRUE} amounts to calling \code{countLexis},
    in which case the arguments \code{new.state}, \code{progressive} and
    \code{precursor.states} will be ignored. }
}
\value{
  A \code{Lexis} object, for which each follow-up interval
  containing the cutpoint is split in two: one before
  and one after the cutpoint. An extra time-scale is added; the time
  since the event at \code{cut}. This is \code{NA} for any follow-up
  prior to the intermediate event.

  The function \code{tsNA20} will replace all missing values in
  timescales with 0. This is commonly meeded when timescales defined as
  time since entry into an inermediate state are used in modeling.
  }
\note{
  The \code{cutLexis} function superficially resembles the
  \code{splitLexis} function. However, the \code{splitLexis} function
  splits on a vector of common cut-points for all rows of the Lexis
  object, whereas the \code{cutLexis} function splits on a single time
  point, which may be distinct for each row, modifies the status
  variables, adds a new timescale and updates the attribute
  "time.since". This attribute is a character vector of the same length
  as the "time.scales" attribute, whose value is '""' if the
  corresponding timescale is defined for any piece of follow-up, and if
  the corresponding time scale is defined by say
  \code{cutLexis(obj,new.state="A",new.scale=TRUE)}, it has the value
  "A".
  }
\details{
  The \code{cutLexis} function allows a number of different ways
  of specifying the cutpoints and of modifying the status variable.

  If the \code{cut} argument is a dataframe it must have columns \code{lex.id},
  \code{cut} and \code{new.state}. The values of \code{lex.id} must be unique.
  In this case it is assumed that each row represents a cutpoint (on the
  timescale indicated in the argument \code{timescale}). This cutpoint will
  be applied to all records in \code{data} with the corresponding \code{lex.id}.
  This makes it possible to apply \code{cutLexis} to a split \code{Lexis} object.

  If a \code{new.state} argument is supplied, the status variable is
  only modified at the time of the cut point. However, it is often
  useful to modify the status variable after the cutpoint when an
  important event occurs. There are three distinct ways of doing this.

  If the \code{progressive=TRUE} argument is given, then a "progressive"
  model is assumed, in which the status can either remain the same or
  increase during follow-up, but never decrease. This assumes that the
  state variables \code{lex.Cst} and \code{lex.Xst} are either numeric or
  ordered factors. In this case, if
  \code{new.state=X}, then any exit status with a value less than
  \code{X} is replaced with \code{X}. The Lexis object
  must already be progressive, so that there are no rows for which the
  exit status is less than the entry status. If \code{lex.Cst} and
  \code{lex.Xst} are factors they must be ordered factors if
  \code{progressive=TRUE} is given.

  As an alternative to the \code{progressive} argument, an explicit
  vector of precursor states, that are considered less severe than the
  new state, may be given. If \code{new.state=X} and
  \code{precursor.states=c(Y,Z)} then any exit status of \code{Y} or
  \code{Z} in the second interval is replaced with \code{X} and all
  other values for the exit status are retained.

  The \code{countLexis} function is a variant of \code{cutLexis} when
  the cutpoint marks a recurrent event, and the status variable is used
  to count the number of events that have occurred. Times given in \code{cut}
  represent times of new events. Splitting with
  \code{countLexis} increases the status variable by 1. If the current
  status is \code{X} and the exit status is \code{Y} before cutting,
  then after cutting the entry status is \code{X}, \code{X+1} for
  the first and second intervals, respectively, and the exit status is
  \code{X+1}, \code{Y+1} respectively. Moreover the values of the status
  is increased by 1 for all intervals for all intervals after the cut
  for the person in question. Hence, a call to \code{countLexis} is
  needed for as many times as the person with most events. But also it
  is immaterial in what order the cutpoints are entered.
}
\author{Bendix Carstensen, Steno Diabetes Center, \email{b@bxc.dk},
        Martyn Plummer, IARC, \email{plummer@iarc.fr}
}
\seealso{
  \code{\link{mcutLexis}},
  \code{\link{addCov.Lexis}},
  \code{\link{splitLexis}},
  \code{\link{Lexis}},
  \code{\link{summary.Lexis}},
  \code{\link{timeSince}},
  \code{\link{boxes.Lexis}}
}
\examples{
# A small artificial example
xx <- Lexis( entry=list(age=c(17,24,33,29),per=c(1920,1933,1930,1929)),
             duration=c(23,57,12,15), exit.status=c(1,2,1,2) )
xx
cut <- c(33,47,29,50)
cutLexis(xx, cut, new.state=3, precursor=1)
cutLexis(xx, cut, new.state=3, precursor=2)
cutLexis(xx, cut, new.state=3, precursor=1:2)
# The same as the last example
cutLexis(xx, cut, new.state=3)

# The same example with a factor status variable
yy <- Lexis(entry = list(age=c(17,24,33,29),per=c(1920,1933,1930,1929)),
            duration = c(23,57,12,15),
            entry.status = factor(rep("alpha",4),
            levels=c("alpha","beta","gamma")),
            exit.status = factor(c("alpha","beta","alpha","beta"),
            levels=c("alpha","beta","gamma")))

cutLexis(yy,c(33,47,29,50),precursor="alpha",new.state="gamma")
cutLexis(yy,c(33,47,29,50),precursor=c("alpha","beta"),new.state="aleph")

## Using a dataframe as cut argument
rl <- data.frame( lex.id=1:3, cut=c(19,53,26), timescale="age", new.state=3 )
rl
cutLexis( xx, rl )
cutLexis( xx, rl, precursor=1 )
cutLexis( xx, rl, precursor=0:2 )

## It is immaterial in what order splitting and cutting is done
xs <- splitLexis( xx, breaks=seq(0,100,10), time.scale="age" )
xs
xsC <- cutLexis(xs, rl, precursor=0 )

xC <- cutLexis( xx, rl, pre=0 )
xC
xCs <- splitLexis( xC, breaks=seq(0,100,10), time.scale="age" )
xCs
str(xCs)
}
\keyword{survival}
