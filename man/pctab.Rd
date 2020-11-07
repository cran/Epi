\name{pctab}
\alias{pctab}
\title{Create percentages in a table}
\description{
Computes percentages and a margin of totals along a given margin of a table.
}
\usage{
pctab(TT, margin = length(dim(TT)), dec=1)
}
\arguments{
  \item{TT}{A table or array object}
  \item{margin}{Which margin should be the the total?}
  \item{dec}{How many decimals should be printed? If 0 or \code{FALSE}
    nothing is printed}
}
\value{
  A table of percentages, where all dimensions except the one specified \code{margin}
  has two extra levels named "All" (where all entries are 100) and "N".
  The function prints the table with \code{dec} decimals.
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://bendixcarstensen.com}.
  }
\seealso{
  \code{\link{addmargins}}
}
\examples{
Aye <- sample( c("Yes","Si","Oui"), 177, replace=TRUE )
Bee <- sample( c("Hum","Buzz"), 177, replace=TRUE )
Sea <- sample( c("White","Black","Red","Dead"), 177, replace=TRUE )
A <- table( Aye, Bee, Sea )
A
ftable( pctab( A ) )
ftable( pctab( addmargins( A, 1 ), 3 ) )
round( ftable( pctab( addmargins( A, 1 ), 3 ), row.vars=3 ), 1)
}
\keyword{ manip }
\keyword{ methods }
\keyword{ array }
