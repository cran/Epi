\name{bootLexis}
\alias{nid}
\alias{nid.Lexis}
\alias{bootLexis}
\title{
Create a bootstrap sample of persons (\code{lex.id}) from a Lexis object 
}
\description{
\code{lex.id} is the identifier of persons in a \code{\link{Lexis}}
object. This is used to sample persons from a Lexis object. If a person
is sampled, all records from this persons is transported to the
bootstrap sample. 
}
\usage{
nid( Lx, \dots )
\method{nid}{Lexis}( Lx, by=NULL, \dots )
bootLexis( Lx, size = NULL, by = NULL )
}
\arguments{
  \item{Lx}{A \code{Lexis} object.}
  \item{\dots}{Parameters passed on to other methods.}
  \item{size}{Numeric. How many persons should be sampled from the
    \code{Lexis} object. Defaults to the number of persons in the
    \code{Lx}, or, if \code{by} is given, to the number of persons in
    each level of \code{by}. If by is given \code{size} can have length
    \code{nid(Lx)}, to indicate how many are sampled from each level of
      \code{by}.} 
  \item{by}{Character. Name of a variable (converted to factor) in the
    \code{Lexis} object.

    Bootstrap sampling is done within each level of by.

    Calculation of the number of persons (\code{lex.id}) is done within
    each level of \code{by}, and a vector returned.
    }
}
\value{
  A Lexis object of the same structure as the input, with \emph{persons}
  bootstrapped. The variable \code{lex.id} has values
  1,2,...,\code{nid(Lx)}.

  \code{nid} counts the number of persons in a Lexis object, possibly by
  \code{by}. If \code{by} is given, a named vector is returned.
}
\author{Bendix Carstensen, \url{BendixCarstensen.com}.}
\seealso{\code{\link{Relevel.Lexis}}}
\examples{
# A small bogus cohort
xcoh <- data.frame( id = c("A", "B", "C"),
                 birth = c("14/07/1952", "01/04/1954", "10/06/1987"),
                 entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                  exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                  fail = c(1, 0, 1),
                   sex = c("M","F","M") )
xcoh <- data.frame( id = c("A", "B", "C"),
                 birth = as.Date(c("14/07/1952","01/04/1954","10/06/1987"),format="\%d/\%m/\%Y"),
                 entry = as.Date(c("04/08/1965","08/09/1972","23/12/1991"),format="\%d/\%m/\%Y"),
                  exit = as.Date(c("27/06/1997","23/05/1995","24/07/1998"),format="\%d/\%m/\%Y"),
                  fail = c(1, 0, 1),
                   sex = c("M","F","M") )
xcoh <- cal.yr( xcoh )
Lcoh <- Lexis( entry = list( per=entry ),
                exit = list( per=exit, age=exit-birth ),
         exit.status = fail,
                data = xcoh )
Lx <- Lcoh[sample(1:3,10,replace=TRUE),]
Lx$per <- Lx$per + runif(10,0,10)
Lx$lex.id <- 1:10
Lx <- splitLexis( Lx, breaks=0:10*10, "age" )
Lx
nid( Lx )
nid( Lx, by="sex" )
bootLexis( Lx )
nid( bootLexis( Lx, size=7 ) )
Li <- bootLexis( Lx, by="id" ) # superfluous
L2 <- bootLexis( Lx, by="sex", size=c(2,5) )
nid( L2, by="sex" )
}
\keyword{manip}
