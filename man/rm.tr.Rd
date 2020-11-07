\name{rm.tr}
\alias{rm.tr}
\title{
Remove transitions from a Lexis object.
}
\description{
Sometimes certain transitions are not of interest. This function removes
these and assigns the risk time in the target state of the transitions
to the originating state.
}
\usage{
rm.tr(obj, from, to)
}
\arguments{
  \item{obj}{
A \code{Lexis} object.
}
  \item{from}{
Character; name of the state from which the transition to be purged
originates. Must be a valid state name for \code{obj}.
}
  \item{to}{
Character; name of the state to which the transition to be purged
targets. Must be a valid state name for \code{obj}.
}
}
\details{
The function removes all transitions from \code{from} to \code{to}, and
assigns all risk time in the \code{to} state after the transition
(\code{lex.dur}) to the \code{from} state. This is only done for risk
time in \code{to} occurring directly after \code{from}. Risk time in
\code{to} occurring after a transition from states different from
\code{from} is not affected. Transitions from \code{to} to another
state, \code{other}, say, will be changed to transitions from
\code{from} to \code{other}. 
}
\value{
A \code{\link{Lexis}} object with the indicated transition removed.
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}.
}
\seealso{
\code{\link{Relevel}}
}
\examples{
data(DMlate)
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit = list(Per=dox),
        exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
               data = DMlate )

# A small subset for illustration
dml <- subset( dml, lex.id \%in\% c(13,15,20,28,40) )

# Cut the follow-up at start of insulin therapy
dmi <- cutLexis( dml, cut = dml$doins,
                      pre = "DM",
                new.state = "Ins" )[,1:10]

# How does it look?
dmi

# Remove all transitions DM -> Ins
rm.tr( dmi, "DM", "Ins" )
}
\keyword{manip}
