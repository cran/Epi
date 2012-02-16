\name{foreign.Lexis}
\Rdversion{1.1}
\alias{msdata}
\alias{msdata.Lexis}
\alias{etm}
\alias{etm.Lexis}
\title{Create a data structures suitable for use with packages mstate, etm.
  }
\description{
The \code{mstate} package requires input in the form of a stacked
dataset with specific variable names. This is provided by
\code{msdata.Lexis}. The resulting dataframe contains the same
information as the result of a call to \code{\link{stack.Lexis}}.

The \code{etm} package requires input (almost) in the form of a
\code{Lexis} object, but with specific column names etc. This is
provided by \code{etm.Lexis}.
}
\usage{
msdata(obj, ...)
\method{msdata}{Lexis}(obj,
                time.scale = timeScales(obj)[1],
                       ... )
\method{etm}{Lexis}( obj,
              time.scale = timeScales(obj)[1],
               cens.name = "cens",
                       s = 0,
                       t = "last",
              covariance = TRUE,
                delta.na = TRUE,
                     ... )
}
\arguments{
  \item{obj}{A \code{\link{Lexis}} object.}
  \item{time.scale}{Name or number of timescale in the \code{Lexis}
    object.}
  \item{cens.name}{Name of the code for censoring used by \code{etm}. It
  is only necessary to change this if one of the states in the
  \code{Lexis} object has name "\code{cens}".}
  \item{s}{Passed on to \code{etm}.}
  \item{t}{Passed on to \code{etm}.}
  \item{covariance}{Passed on to \code{etm}.}
  \item{delta.na}{Passed on to \code{etm}.}
  \item{\dots}{Further arguments.}
}
\value{
  \code{msdata.Lexis} returns a dataframe with the \code{Lexis} specific
  variables stripped, and with the following added: \code{id},
  \code{Tstart}, \code{Tstop}, \code{from}, \code{to}, \code{trans},
  \code{status}, which are used in the \code{mstate} package.

  \code{etm.Lexis} transforms the \code{Lexis} object into a dataframe
  suitable for analysis by the function \code{etm} from the \code{etm}
  package, and actually calls this function, so returns an object of
  class \code{etm}.
}
\author{
Bendix Carstensen, \email{bxc@steno.dk}, \url{www.biostat.ku.dk/~bxc}
}
\examples{
data(DMlate)
str(DMlate)
dml <- Lexis( entry=list(Per=dodm,Age=dodm-dobth,DMdur=0),
               exit=list(Per=dox),
        exit.status=factor(!is.na(dodth),labels=c("DM","Dead")),
               data=DMlate )
dmi <- cutLexis( dml, cut=dml$doins, new.state="Ins", pre="DM" )
summary( dmi )
ms.dmi <- msdata.Lexis( dmi )
summary( dmi )
# Check that all the transitions and person-years got across.
with( ms.dmi, rbind( table(status,trans),
                     tapply(Tstop-Tstart,trans,sum) ) )

# Use the etm package directly with a Lexis object
dmi <- subset(dmi,lex.id<1000)
etm.A <- etm.Lexis( dmi, time.scale=2 )
etm.D <- etm.Lexis( dmi, time.scale=3 )
plot( etm.A, col=rainbow(5), lwd=2, lty=1, xlab="Age" )
plot( etm.D, col=rainbow(5), lwd=2, lty=1, xlab="DM duration" )
}
\seealso{
  \code{\link{stack.Lexis}},
  \code{\link[mstate:msprep]{msdata}},
  \code{\link[etm:etm]{etm}} }
\keyword{survival}
