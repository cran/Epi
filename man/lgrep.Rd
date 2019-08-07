\name{lgrep}
\alias{fgrep}
\alias{ngrep}
\alias{lgrep}
\title{
Convenience versions of grep
}
\description{
Often you want the elements of a vector (or its names or levels) that
meet a certain pattern. But \code{grep} only gives you the position, so
these functions are designed to give you that. 
}
\usage{
fgrep( pattern, x, ... )
ngrep( pattern, x, ... )
lgrep( pattern, x, ... )
}
\arguments{
  \item{pattern}{Pattern searched for.}
  \item{x}{Object where \code{pattern} is searched. Or in whose \code{names}
    or \code{levels} attributes \code{pattern} is sought.}
  \item{...}{Arguments passed on to \code{\link[base]{grep}}.}
}
\value{Elements of the input \code{x} (\code{fgrep}) or its names
  attribute (\code{ngrep}) or levels attribute (\code{lgrep}).
}
\author{Bendix Carstensen, \email{b@bxc.dk}, \url{http://BendixCarstensen.com} }
\seealso{\code{\link{grep}}}
\examples{
ff <- factor( ll <- paste( sample( letters[1:3], 20, replace=TRUE ),
                           sample( letters[1:3], 20, replace=TRUE ), sep="" ) )
ff
fgrep( "a", ff )
fgrep( "a", ll )
ngrep( "a", ff )
lgrep( "a", ff )
lgrep( "a", ff, invert=TRUE )
}
\keyword{ manip }
