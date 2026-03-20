\name{LexisScarf}
\alias{LexisScarf}
\title{
Generate a knitting pattern for a tube-shaped scarf with Lexis lines
}
\description{
The function generates a dotted pattern that refers to stitches of a
knitting pattern. Since there are several colors involved the back side
of the knit becomes quite messy, hence the pattern is for a tube-shaped
scarf, that will have all the messiness on the inside. So use a circular
needle for the work. And send the author a picture of your finished scarf.
}
\usage{
LexisScarf(nw = 4,   # fields around (width of scarf)           
           nl = 20,  # fields along  (length)           
           nm = 23,  # no. stitches per field   
           np = 30,  # no. persons              
        scale = 1.5, # maximal FU in units of nw
          clr = c("red", "blue", "limegreen"))
}
\arguments{
  \item{nw}{number of fields around the tube}
  \item{nl}{number of fields along the tube}
  \item{nm}{no. stitches per field}
  \item{np}{no. lifelines (persons)}
  \item{scale}{maximal length of life lines in units of \code{nw}}
  \item{clr}{vector of colors to be used for life lines}
}
\details{The function produces a graph; ask for a plotting area with
  aspect ratio of \code{nl} to \code{nw}.
}
\value{NULL. The function is used for its side-effect of
  producing a graph of the knitting pattern.
}
\author{Bendix Carstensen,
  \url{http://bendixcarstensen.com},
  \email{b@bxc.dk}}
\examples{
# pdf("scarf.pdf", height = 10, width = 50)
LexisScarf()
# dev.off()
}
\keyword{manip}
