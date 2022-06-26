\name{Lexis2msm}
\alias{Lexis2msm}
\title{Convert a Lexis obejct to a data set suitable for input to the
  \code{msm:msm} function.} 
\description{The number of records in the resulting dataset will have a
  number of records that is normally \code{nrec(Lx) + \link{nid}(Lx)}, that is
  one extra record for each person. If there are 'holes' in persons'
  follow-up, each hole will also generate an extra record in the result. 
}
\usage{
Lexis2msm(Lx,
       state = "state",
     verbose = FALSE)
}
\arguments{
  \item{Lx}{
A \code{\link{Lexis}} object.
}
  \item{state}{
Character; the name of the state variable in the result.
}
  \item{verbose}{
If true, you will be reminded what the function did.
}
}
\value{
A data frame of class \code{msmLexis} with the timescales preserved and \code{lex.id}
preserved but with other \code{lex.} variables removed. 

Has more records than the original \code{Lexis} object
}
\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link{Lexis}}
}
\examples{
example(mcutLexis)
# we now have the Lexis object L3:
summary(L3)

# data frame for use with msm
msm3 <- Lexis2msm(L3)

# see the difference 
subset(  L3, lex.id \%in\% 1:3)
subset(msm3, lex.id \%in\% 1:3)
timeScales(msm3)
}
\keyword{survival}
\keyword{manip}
