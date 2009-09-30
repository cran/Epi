\name{plotevent}
\alias{plotevent}
\title{ Plot Equivalence Classes }
\description{
  For interval censored data, segments of times between last.well and first.ill are plotted for each conversion in the data. It also plots the equivalence classes.
}
\usage{
plotevent(last.well, first.ill, data)
}
\arguments{
  \item{last.well}{ Time at which the individuals are
    last seen negative for the event }
  \item{first.ill}{ Time at which the individuals are
    first seen positive for the event }
  \item{data}{ Data with a transversal shape }
}
\details{
  last.well and first.ill should be written as character in the function.
}
\value{
  Graph
}
\references{ 
Carstensen B. Regression models for interval censored survival data:
application to HIV infection in Danish homosexual men.Stat Med. 1996 Oct
30;15(20):2177-89. 

Lindsey JC, Ryan LM. Tutorial in biostatistics methods for
interval-censored data.Stat Med. 1998 Jan 30;17(2):219-38. 
  }
\author{
  Delphine Maucort-Boulch, Bendix Carstensen, Martyn Plummer
  }
\seealso{
  \code{\link{Icens}}
  }
% \examples{
% }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
