\name{legendbox}
\alias{legendbox}
\title{Draw a box with text explaining the numbers in and between boxes
  from \code{boxes.MS} and  \code{boxes.Lexis}
}
\description{When drawing boxes describing a multistate model a legend
  explaining the numbers in the plot is required. \code{legendbox} does
  this. 
}
\usage{
legendbox(x, y,
         state = "State",
            py = "Person-time",
         begin = "no. begin",
           end = "no. end",
         trans = "Transitions",
         rates = "\n(Rate)",
          font = 1,
         right = !left,
          left = !right,
           ...)
}
\arguments{
  \item{x}{
x-coordinate of the center of the box.
}
  \item{y}{
y-coordinate of the center of the box.
}
  \item{state}{
Text describing the state
}
  \item{py}{
Text describing the risk time
}
  \item{begin}{
Text describing the no. persons starting FU in state
}
  \item{end}{
Text describing the no. persons ending FU in state
}
  \item{trans}{
Text describing the no. of transitions
}
  \item{rates}{
Text describing the rates
}
  \item{font}{
Font to use for the text 
}
  \item{right}{
Should a text describing arrow texts be on the r.h.s. of the box?
Defaults to TRUE.
}
  \item{left}{
Should a text describing arrow texts be on the l.h.s. of the box?
}
  \item{\dots}{
Arguments passed on to \code{tbox}
}
}
\value{None.
}
\details{The function is called for its side effect of adding an
  explanatory box to the plot. If \code{right} is true, an explanation
  of events and rates are added to the right of the box. Similarly for
  \code{left}. It is admissible that \code{left == right}.
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link{boxes.Lexis}}
}
\keyword{aplot}
