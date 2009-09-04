\name{mstate.Lexis}
\Rdversion{1.1}
\alias{mstate.Lexis}
\title{Create a dataframe suitable for use with the mstate package.
  }
\description{
The \code{mstate} package requires input in the form of a stacked
dataset with specific variable names. This is provided by this
function. The resulting dataframe contains the same information as the
result of a call to \code{\link{stack.Lexis}}.
}
\usage{
mstate.Lexis(obj, time.scale = timeScales(obj)[1])
}
\arguments{
  \item{obj}{A \code{\link{Lexis}} object.}
  \item{time.scale}{Name or number of timescale in the \code{Lexis} object.}
}
\value{
  A dataframe with the \code{Lexis} specific variables stripped, and
  with the following added: \code{id}, \code{Tstart}, \code{Tstop},
  \code{from}, \code{to}, \code{trans}, \code{status}, which are used in
  the function \code{\link[mstate]{mstate}} from the \code{mstate} package.
}
\author{
Bendix Carstensen, \email{bxc@steno.dk}, \url{www.biostat.ku.dk/~bxc}
}
\examples{
data(DMlate)
str(DMlate)
dml <- Lexis( entry=list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
               exit=list(Per=dox),
        exit.status=factor(!is.na(dodth),labels=c("DM","Dead")),
               data=DMlate )
dmi <- cutLexis( dml, cut=dml$doins, new.state="Ins", pre="DM" )
summary( dmi )
ms.dmi <- mstate.Lexis( dmi )
summary( dmi )
# Check that all the transitions and person-years got across.
with( ms.dmi, rbind( table(status,trans),
                    tapply(Tstop-Tstart,trans,sum) ) )
}
\seealso{\code{\link{stack.Lexis}}}
\keyword{survival}
