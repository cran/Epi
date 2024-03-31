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
          font = 1, ...)
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
  \item{\dots}{
Argumants passed on to \code{tbox}
}
}
\value{None.
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link{boxes.Lexis}}
}
\keyword{aplot}
