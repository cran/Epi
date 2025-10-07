\name{xcutLexis}
\alias{xcutLexis}
\title{
Cut follow-up in a \code{\link{Lexis}} object by event date(s) while
preserving the original states. This is essentially a
cross-classification of the original states and the new ones, hence the "\code{x}".
}
\description{
With a multistate \code{\link{Lexis}} object we might want a
classification of the follow-up according to a set of events independent
of the states in the Lexis object, while keeping the the original states.
This is what \code{xcutLexis}
}
\usage{
xcutLexis(Lx, cut, timescale = 1, sep = ".")
}
\arguments{
  \item{Lx}{
A \code{Lexis} object
}
  \item{cut}{
    A data frame with at most one row per \code{lex.id}, columns must be:
    \itemize{
    \item \code{lex.id}: person-id referring to \code{lex.id} in \code{Lx}
    \item \code{cut}, a numeric; times of event om the time scale \code{timescale}
    \item \code{new.state}, character or factor; name of the new state
    occurring at time \code{cut} 
            }
}
  \item{timescale}{
The timescale that \code{cut} refers to. Numeric or character.
}
  \item{sep}{
Character; string used in forming new state names by joining original
state names in \code{Lx} and state names in \code{cut$new.states} 
}
}
\details{
The function was motivated by a follow-up through states of
multimorbidity with state names \code{0morb},  \code{1morb},
\code{2morb}, etc. where we wanted a subdivision of each state by the
presence of T2 diabetes, resulting in states
\code{0morb},  \code{1morb}, \code{2morb},
\code{0morb.T2},  \code{1morb.T2}, \code{2morb.T2}.

At most 1 transition per person is allowed in \code{cut}.
}
\value{
A \code{Lexis} object with states as in \code{Lx} plus states named by
concatenating \code{Lx} state names with names in \code{cut$new.state}
}
\author{Bendix Carstensen, Steno Diabetes Center Copenhagen, \email{b@bxc.dk},
\url{http://bendixcarstensen.com}
}
\seealso{
  \code{\link{mcutLexis}},
  \code{\link{rcutLexis}},
  \code{\link{addCov.Lexis}},
  \code{\link{splitLexis}},
  \code{\link{Lexis}},
  \code{\link{summary.Lexis}},
  \code{\link{timeSince}},
  \code{\link{boxes.Lexis}}
}
\examples{
example(DMlate)
levels(dmi)
#
# show transitions between states in dmi
boxes(dmi, boxpos = list(x = rep(50, 4),
                         y = c(65, 35, 95, 5)),
           scale.R = 1000,
           show.BE = TRUE)
#
# randomly generated intermediate events X and Y
set.seed(1952)
cutXY <- data.frame(lex.id = unique(dmi$lex.id),          # one row per id
                       cut = runif(nid(dmi), 1995, 2008), # event dates
                 new.state = sample(c("X","Y"),           # event types
                                    nid(dmi),             # only 4 in 10 has an event
                                    repl = TRUE))[runif(nid(dmi)) < 0.4, ]
#
# cut at these event dates but also keep original states
Lxy <- xcutLexis(dmi, cutXY)
levels(Lxy)
#
# reorder the levels for easier specification of box places
Lxy <- Relevel(Lxy, as.vector(t(outer(levels(dmi),
                                      c("", ".X", ".Y"),
                                      paste0))))
levels(Lxy)
summary(Lxy)
#
# resulting transitions between boxes
boxes(Lxy, boxpos = list(x = rep(c(50, 15, 85), 4),
                         y = rep(c(65,35,95,5), each = 3)),
           scale.R = 1000,
           show.BE = TRUE,
           cex = 0.8)
}
\keyword{survival}
\concept{Lexis}
\concept{multistate models}
\concept{transition rates}
\concept{state transition}
