\name{Relevel}
\alias{Relevel}
\alias{Relevel.factor}
\title{Reorder and combine levels of a factor}
\description{
  The levels of a factor are re-ordered so that the levels specified by
  \code{ref} appear first and remaining levels are moved down. This is
  useful for \code{contr.treatment} contrasts which take the first level
  as the reference. Factor levels may also be combined; two possibilities for
  specifying this are supported: hard coding or table look-up.
}
\usage{
\method{Relevel}{factor}( x, ref, first = TRUE, collapse="+",
                 xlevels=TRUE, nogroup=TRUE, \dots )
}
\arguments{
  \item{x}{A(n unordered) factor}
  \item{ref}{Vector, list or data frame, array, matrix or table.

    If \code{ref} is a vector (integer or character), it is assumed it
    contains the names or numbers of levels to be the first ones; non
    mentioned levels are kept.

    If \code{ref} is a list (but not a data frame), factor levels
    mentioned in each list element are combined. If the list is named
    the names are used as new factor levels, otherwise new level names
    are constructed from the old.

    If \code{ref} is a data frame or 2-dimensional array, matrix or
    table, the first column is assumed to have unique levels of \code{x}
    and the second to have groupings of this, respectively.

    }
  \item{first}{Should the levels mentioned in \code{ref} (if it is a
    list) come before those not?} 
  \item{collapse}{String used when constructing names for combined factor levels.}
  \item{xlevels}{Logical. Should all levels in the 2nd column of
    \code{ref} be maintained as levels of the result, or (if
    \code{FALSE}) only the actually occurring.}
  \item{nogroup}{Logical. Should levels present in the input but not in
    the 1st column of \code{ref} be maintained as levels after the
    grouping? If \code{FALSE}), NAs will be returned for such elements.}
  \item{\dots}{Arguments passed on to other methods.}
}
\value{
  An unordered factor, where levels of \code{x} have been reordered
  and/or collapsed.
}
\details{
  The facility where \code{ref} is a two-column matrix mimics the
  SAS-facility of formats where a dataset can be used to construct a
  format --- SAS format is the grouping tool for variable values. Using
  this facility will preserve the order of levels if \code{ref[,2]} is a
  factor.
}
\author{Bendix Carstensen, \url{BendixCarstensen.com}, Lars Jorge Diaz}
\seealso{\code{\link{Relevel.Lexis}}}
\examples{
# Grouping using a list (hard coding)
#
ff <- factor( sample( letters[1:5], 100, replace=TRUE ) )
table( ff, Relevel( ff, list( AB=1:2, "Dee"=4, c(3,5) ) ) )
table( ff, Relevel( ff, list( 5:4, Z=c("c","a") ), coll="-und-",first=FALSE ) )

## Grouping using a two-column matrix as input:
## A factor with levels to be grouped together
ff <- factor( c("Bear","Bear","Crocodile","Snake","Crocodile","Bear") )
ff
## A grouping table
( gg <- data.frame( Animal=c("Bear","Whale","Crocodile","Snake","Eagle"),
                     Class=c("Mammal","Mammal","Reptile","Reptile","Bird") ) )
Relevel( ff, gg, xlevels=FALSE )
Relevel( ff, gg )
Relevel( ff, gg[c(1:5,5:1),] )

## This crashes with an error
try( Relevel( ff, rbind( gg, c("Bear","Reptile") ) ) )

ff <- factor(c(as.character(ff),"Jellyfish","Spider"))
Relevel( ff, gg )
Relevel( ff, gg, xlevels=FALSE )
# If you do not want unknown animals classified, this returns NAs:
Relevel( ff, gg, nogroup=FALSE )
# Will give error because re-classification is ambiguous
( GG <- rbind(gg,c("Eagle","Mammal"),c("Snake","Mammal")) )
try( Relevel( ff, GG ) )
}
\keyword{manip}
