\name{Relevel}
\alias{Relevel}
\alias{Relevel.factor}
\title{Reorder and combine levels of a factor}
\description{
  The levels of a factor are re-ordered so that the levels specified by
  \code{ref} is first and the others are moved down. This is useful for
  \code{contr.treatment} contrasts which take the first level as the
  reference. Levels may also be combined.
}
\usage{
\method{Relevel}{factor}( x, ref, first = TRUE, collapse="+", \dots )
}
\arguments{
  \item{x}{An unordered factor}
  \item{ref}{The names or numbers of levels to be the first. If \code{ref}
  is a list, factor levels mentioned in each list element are
  combined. If the list is named the names are used as new factor levels.}
  \item{first}{Should the levels mentioned in ref come before those not?}
  \item{collapse}{String used when collapsing factor levels.}
  \item{\dots}{Arguments passed on to other methods.}
}
\value{
  An unordered factor, where levels of \code{x} have been reordered
  and/or collapsed.
}
\seealso{\code{\link{Relevel.Lexis}}}
\examples{
ff <- factor( sample( letters[1:5], 100, replace=TRUE ) )
table( ff, Relevel( ff, list( AB=1:2, "Dee"=4, c(3,5) ) ) )
table( ff, rr=Relevel( ff, list( 5:4, Z=c("c","a") ), coll="-und-", first=FALSE ) )
}
\keyword{manip}
