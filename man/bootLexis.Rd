\name{bootLexis}
\alias{nid}
\alias{nid.Lexis}
\alias{bootLexis}
\title{
Create a bootstrap sample of persons (as identified by \code{lex.id}) from a Lexis object 
}
\description{
\code{lex.id} is the person identifier in a \code{\link{Lexis}}
object. This is used to sample persons from a Lexis object. If a person
is sampled, all records from this persons is transported to the
bootstrap sample. 
}
\usage{
nid( Lx, \dots )
\method{nid}{Lexis}( Lx, by=NULL, \dots )
bootLexis( Lx, size = NULL, by = NULL, replace=TRUE )
}
\arguments{
  \item{Lx}{A \code{Lexis} object.}
  \item{\dots}{Parameters passed on to other methods.}
  \item{size}{Numeric. How many persons should be sampled from the
    \code{Lexis} object. Defaults to the number of persons in the
    \code{Lx}, or, if \code{by} is given, to the number of persons in
    each level of \code{by}. If \code{by} is given, \code{size} can have length
    \code{length(unique(by))}, to indicate how many are sampled from
    each level of \code{by}.}  
  \item{by}{Character. Name of a variable (converted to factor) in the
    \code{Lexis} object.

    Bootstrap sampling is done within each level of by.

    Calculation of the number of persons (\code{lex.id}) is done within
    each level of \code{by}, and a vector returned.
  }
  \item{replace}{Should persons be sampled by replacement? Default is
    \code{TRUE}. Setting \code{replace} to \code{FALSE} enables
    selecting a random subset of persons from the Lexis object.}
}
\value{  
  \code{bootLexis} returns a Lexis object of the same structure as the
    input, with \emph{persons} bootstrapped. The variable \code{lex.id}
    in the resulting \code{Lexis} object has values 1,2,... The original
    values of \code{lex.id} from \code{Lx} are stored in the variable
    \code{old.id}.

  \code{nid} counts the number of persons in a Lexis object, possibly by
  \code{by}. If \code{by} is given, a named vector is returned.
}
\author{Bendix Carstensen, \url{http://bendixcarstensen.com}.}
\seealso{\code{\link{Relevel.Lexis}},\code{\link{subset.Lexis}}}
\examples{
# A small bogus cohort
xcoh <- data.frame( id = c("A", "B", "C"),
                 birth = c("1952-07-14", "1954-04-01", "1987-06-10"),
                 entry = c("1965-08-04", "1972-09-08", "1991-12-23"),
                  exit = c("1997-06-27", "1995-05-23", "1998-07-24"),
                  fail = c(1, 0, 1),
                   sex = c("M","F","M") )
# Convert to calendar years
for( i in 2:4 ) xcoh[,i] <- cal.yr(xcoh[,i])
xcoh
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
Lb <- bootLexis( Lx )
head( Lb )
nid( bootLexis( Lx, size=7 ) )
Li <- bootLexis( Lx, by="id" ) # superfluous
summary( Lx )
summary( Li )
L2 <- bootLexis( Lx, by="sex", size=c(2,5) )
nid( L2, by="sex" )
summary( L2 )
}
\keyword{manip}
