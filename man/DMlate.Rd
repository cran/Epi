\name{DMlate}
\Rdversion{1.1}
\alias{DMlate}
\alias{DMrand}
\docType{data}
\title{
The Danish National Diabetes Register.
}
\description{
These two datasets each contain a random sample of 10,000 persons from
the Danish National Diabetes Register. \code{DMrand} is a random sample
from the register, whereas \code{DMlate} is a random sample among those
with date of diagnosis after 1.1.1995.
}
\usage{data(DMrand)
       data(DMlate)}
\format{
  A data frame with 10000 observations on the following 6 variables.
  \describe{
    \item{\code{sex}}{Sex, a factor with levels \code{M} \code{F}}
    \item{\code{dobth}}{Date of birth}
    \item{\code{dodm}}{Date of inclusion in the register}
    \item{\code{dodth}}{Date of death}
    \item{\code{doins}}{Date of first insulin prescription}
    \item{\code{dox}}{Date of exit from follow-up.}
  }
}
\details{All dates are given in fractions of years, so 1997.00
  corresponds to 1 January 1997 and 1997.997 to 31 December 1997.
}
\source{
Danish National Board of Health.
}
\references{
B Carstensen, JK Kristensen, P Ottosen and K Borch-Johnsen:
The Danish National Diabetes Register: Trends in incidence, prevalence and
mortality, Diabetologia, 51, pp 2187--2196, 2008.

In partucular see the appendix at the end of the paper.
}
\examples{
data(DMlate)
str(DMlate)
dml <- Lexis( entry=list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit=list(Per=dox),
        exit.status=factor(!is.na(dodth),labels=c("DM","Dead")),
               data=DMlate )
# Split follow-up at Insulin
dmi <- cutLexis( dml, cut=dml$doins, new.state="Ins", pre="DM" )
summary( dmi )
# Introduce a new timescale
dmi <- cutLexis( dml, cut=dml$doins, new.state="Ins", pre="DM", new.scale=TRUE )
head( dmi )
# Split the states following insulin and explictily name the new timescale
dmi <- cutLexis( dml, cut=dml$doins, new.state="Ins",
                      pre="DM", new.scale="Instime", split.states=TRUE )
summary( dmi )
}
\keyword{datasets}
