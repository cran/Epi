\name{rcutLexis}
\alias{rcutLexis}
\title{
A function to cut follow-up at intermediate event times.
}
\description{
Cuts follow-up at intermediate event times, multiple events per person
are allowed, as well as recurrences of the sme type of event. The
resulting states only refer to the last assumed state, unlike the result
from \code{\link{mcutLexis}}.
}
\usage{
rcutLexis( Lx, cut,
    timescale = 1,
    precursor.states = transient(Lx))
}
\arguments{
  \item{Lx}{
  A \code{\link{Lexis}} object to be amended,.
}
  \item{cut}{
  A data frame with columns \code{lex.id}, \code{cut} (event times) and
  \code{new.state} (event type, character)
}
  \item{timescale}{
  What time scale do values in \code{cut$cut} refer to. Numeric or character.
}
  \item{precursor.states}{an optional vector of states to be considered
    as "less severe" than \code{new.state}. See Details in the
    documentation of \code{\link{cutLexis}}
}
}
\value{A \code{\link{Lexis}} object with follow-up cut at the event
  times supplied in \code{cut}}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
  \code{\link{cutLexis}},
  \code{\link{mcutLexis}},
  \code{\link{addCov.Lexis}},
  \code{\link{Lexis}},
  \code{\link{splitLexis}}
}
\examples{
df <- data.frame(lex.id = rep(c(3, 7), c(3, 5)))
df$new.state <- sample(LETTERS[2:4], 8, r = TRUE) 
df$cut <- round(runif(8) * 100) + 1
df

Lx <- Lexis( exit = list(time=c(89, 97)),
               id = c(3, 7),
      exit.status = factor(c("A", "X")) )
Lx 
rcutLexis(Lx, df, pre = "A")
        }
\keyword{survival}
