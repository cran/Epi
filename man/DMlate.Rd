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
with date of diagnosis after 1.1.1995. All dates are radomly jittered by
adding a U(-7,7) (days).
}
\usage{data(DMrand)
       data(DMlate)}
\format{
  A data frame with 10000 observations on the following 7 variables.
  \describe{
    \item{\code{sex}}{Sex, a factor with levels \code{M} \code{F}}
    \item{\code{dobth}}{Date of birth}
    \item{\code{dodm}}{Date of inclusion in the register}
    \item{\code{dodth}}{Date of death}
    \item{\code{dooad}}{Date of 2nd prescription of OAD}
    \item{\code{doins}}{Date of 2nd insulin prescription}
    \item{\code{dox}}{Date of exit from follow-up.}
  }
}
\details{All dates are given in fractions of years, so 1998.000
  corresponds to 1 January 1998 and 1998.997 to 31 December 1998.

  All dates are randomly perturbed by a small amount, so no real
  persons have any of the combinations of the 6 dates in the
  dataset. But results derived from the data will be quite close to
  those that would be obtained if the entire 'real' diabetes register
  were used.
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

# Cut the follow-up at insulin start, and introduce a new timescale,
# and split non-precursor states
system.time(
dmi <- cutLexis( dml, cut = dml$doins,
                      pre = "DM",
                new.state = "Ins",
                new.scale = "t.Ins",
             split.states = TRUE ) )
summary( dmi )
}
\keyword{datasets}
