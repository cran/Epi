\name{mat2pol}
\alias{mat2pol}
\title{Plot columns of a matrix as stacked areas.
}
\description{\code{mat}rix to \code{pol}ygon: Plot columns of a matrix
  as stacked areas.
}
\usage{
mat2pol( pm,
       perm = 1:ncol(pm),
          x = as.numeric(rownames(pm)),
        col = rainbow(ncol(pm)),
         yl = 0:1,
     append = FALSE,
        ... )
}
\arguments{
  \item{pm}{Numerical matrix.}
  \item{perm}{integer vector of length \code{ncol(pm)}, used to permute
    the columns of \code{pm}.}
  \item{x}{Numeric. The x-axis of the plot.}
  \item{col}{Colors of the areas.}
  \item{yl}{y-axis limits.}
  \item{append}{Logical. Should the polygons be added to an exiating plot}
  \item{\dots}{Further parameters passed to \code{plot}.}
}
\details{The function is originally intended to plot stacked
  probabilities, hence the default of \code{0:1} for the y-axis.
}
\value{A matrix of \code{ncol(pm)+1} columns with the first equal to 0,
  and the remaining the cumulative sum of the columns of
  \code{pm[perm]}.  

  The function is called for its side effect - the stacked polygons.
}
\author{Bendix Carstensen
}
\examples{
M <- cbind( sort(runif(10)), sort(runif(10)), sort(runif(10)) )
pm <- sweep( M, 1, apply(M,1,sum), "/" )
mat2pol( pm )
}
\keyword{manip}
