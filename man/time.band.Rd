\name{timeBand}
\alias{timeBand}
\alias{breaks}
\title{Extract time band data from a split Lexis object}
\description{
  The break points of a \code{Lexis} object (created by a call to
  \code{splitLexis}) divide the follow-up intervals into time bands
  along a given time scale. The \code{breaks} function returns
  the break points, for a given time scale, and the \code{timeBand}
  classifies each row (=follow-up interval) into one of the time bands.
}
\usage{
timeBand(lex, time.scale, type="integer")
breaks(lex, time.scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lex}{an object of class \code{Lexis}}
  \item{time.scale}{a character or integer vector of length 1
    identifying the time scale of interest}
  \item{type}{a string that determines how the time bands are labelled.
    See Details below}
}
\details{
  Time bands may be labelled in various ways according to the
  \code{type} argument. The permitted values of the \code{type}
  argument, and the corresponding return values are:
  \describe{
  \item{"integer"}{a numeric vector with integer codes starting from 0.}
  \item{"factor"}{a factor (unordered) with labels "(left,right]"}
  \item{"left"}{the left-hand limit of the time band}
  \item{"middle"}{the midpoint of the time band}
  \item{"right"}{the right-hand limit of the time band}
  }
}    
\value{
  The \code{breaks} function returns a vector of break points
  for the \code{Lexis} object, or NULL if no break points have been
  defined by a call to \code{splitLexis}.  The \code{timeBand}
  function returns a numeric vector or factor, depending on the value
  of the \code{type} argument.
}
\author{Martyn Plummer}
\note{
  A newly created \code{Lexis} object has no break points defined.
  In this case, \code{breaks} will return NULL, and
  \code{timeBand} will a vector of zeros.
}  
\examples{
data(diet)
diet <- cal.yr(diet)
diet.lex <- Lexis(entry=list(period=doe),
                   exit=list(period=dox, age=dox-dob),
            exit.status=chd,
                   data=diet)
diet.split <- splitLexis(diet.lex, breaks=seq(40,70,5), "age" )
age.left <- timeBand(diet.split, "age", "left")
table(age.left)
age.fact <- timeBand(diet.split, "age", "factor")
table(age.fact)
age.mid <- timeBand(diet.split, "age", "mid")
table(age.mid)
}
\seealso{\code{\link{Lexis}}}
\keyword{attribute}
