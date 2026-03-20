\name{overlap.Lexis}
\alias{overlap}
\alias{overlap.Lexis}
\title{Check if records in a Lexis object represent overlapping epochs}
\description{
When constructing a \code{\link{Lexis}} object from a data frame with
multiple records there is no guarantee that the records from a given
person will represent disjoint time intervals (epochs), \code{overlap}
checks this.
}
\usage{
\method{overlap}{Lexis}(Lx, \dots)
}
\arguments{
  \item{Lx}{A \code{Lexis} object}
  \item{\dots}{ignored}
}
\value{
If no overlapping epochs are found, a message is printed and \code{NULL}
is returned. Otherwise a vector of unique ids (from \code{Lx$lex.id})
for the persons with overlapping records is returned invisibly
}
\author{Bendix Carstensen,
  \url{http://bendixcarstensen.com},
  \email{b@bxc.dk}}
\seealso{\code{\link{Lexis}}}
\examples{
# run the DMlate example to get the Lexis object dmi with
# some persons with multiple (2 that is) records:
example(DMlate)
overlap.Lexis(dmi)
#
# now find 4 persons with multiple records
tt <- table(dmi$lex.id)
has2 <- dmi$lex.id %in% names(tt[tt > 1])[1:4]
table(has2)
#
# add time (0.3 y) to lex.dur for these persons to produce overlapping intervals
dmi$lex.dur[has2] <- dmi$lex.dur[has2] + 0.3
overlap.Lexis(dmi)
}
\keyword{manip}
\concept{Lexis checks}
