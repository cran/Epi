\name{NArray}
\alias{NArray}
\alias{ZArray}
\title{Set up an array of NAs, solely from the list of dimnames
}
\description{Defines an array of NAs, solely from the list of dimnames
}
\usage{
NArray( x, cells=NA )
ZArray( x, cells=0 )
}
\arguments{
  \item{x}{A (possibly named) list to be used as dimnames for the
    resulting array}
  \item{cells}{Value(s) to fill the array}
}
\details{This is a simple useful way of defining arrays to be used for
  collection of results. The point is that everything is defined from
  the named list, so in the process of defining what you want to
  collect, there is only one place in the program to edit. It's just a
  wrapper for \code{array}. \code{ZArray} is just a wrapper for
  \code{NArray} with a different default.
}
\value{An array with \code{dimnames} attribute \code{x}, and all values
  equal to \code{cells}.
}
\author{Bendix Carstensen
}
\examples{
ftable(
NArray( list(Aye = c("Yes", "Si", "Oui"),
             Bee = c("Hum", "Buzz"),
             Sea = c("White", "Black", "Red", "Dead") ) ) )
}