\name{AaJ.Lexis}
\alias{AaJ.Lexis}
\title{The Aalen-Johansen estimator of state probabilities from a
  multistate \code{Lexis} object.
}
\description{
The Aalen-Johansen estimator is computed on the basis of a
\code{\link{Lexis}} multistate object along a given time scale. The
function is merely a wrapper for the \code{\link[survival]{survfit}}.
}
\usage{
\method{AaJ}{Lexis}(Lx,
               formula = ~ 1,
             timeScale = 1)
}
\arguments{
  \item{Lx}{A \code{Lexis} object. The starting state must be the first
    among \code{levels(Lx)}.
}
  \item{formula}{A one-sided formula passed on to \code{survfit}.
}
  \item{timeScale}{Chararter or integer, selecting one of the timescales
    of the \code{Lexis} object.
}
}
\value{An object of class \code{survfitms} --- see
  \code{\link[survival]{survfit}}. 
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link[survival]{survfit}}
\code{\link{ci.Crisk}}
}
\examples{
data(DMlate)
str(DMlate)
dml <- Lexis(entry = list(Per = dodm,
                          Age = dodm-dobth,
                        DMdur = 0 ),
              exit = list(Per = dox),
       exit.status = factor(!is.na(dodth),
                            labels = c("DM","Dead")),
              data = DMlate )

# Cut the follow-up at insulin start
dmi <- cutLexis(dml,
                cut = dml$doins,
          new.state = "Ins",
        split.state = TRUE)
summary( dmi )
ms <- AaJ.Lexis(dmi, timeScale = "DMdur")
class(ms)
ms$states
head(ms$pstate)
}
\keyword{models}
