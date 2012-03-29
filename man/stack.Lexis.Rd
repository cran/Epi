\name{stack.Lexis}
\Rdversion{1.1}
\alias{stack.Lexis}
\alias{tmat}
\alias{tmat.Lexis}
\title{
Functions to facilitate analysis of multistate models.
}
\description{
\code{stack.Lexis} produces a stacked object suited for analysis of
several transitions simultaneously.
}
\usage{
\method{stack}{Lexis}(x, ...)
tmat( x, ... )
\method{tmat}{Lexis}(x, Y=FALSE, mode = "numeric", ...)
}
\arguments{
  \item{x}{A \code{\link{Lexis}} object.}
  \item{Y}{Logical. Should the risk time be put in the diagonal? This is
    a facility which is used by \code{\link{boxes.Lexis}}.}
  \item{mode}{Should the matrix be returned as a numeric matrix with
    \code{NA}s at unused places or (\code{mode="logical"}) as a logical
    matrix with \code{FALSE} on the diagonal.}
  \item{\dots}{Not used.}
}
\value{
\code{tmat.Lexis} returns a square transition matrix, classified by the
levels of \code{lex.Cst} and \code{lex.Xst}, for every transition
occurring the entry is the number of transitions occurring and \code{NA}
in all oter entries. If \code{Y=TRUE}, the diagonal will contain the
risk time in each of the states.

\code{stack.Lexis} returns a dataframe to be used for analysis of
multistate data when all transitions are modelled together, for example
if some parameters are required to be the same for different transitions.

The dataframe has same variables as the original \code{Lexis} object,
but with each record duplicated as many times as there are possible
exits from the current state, \code{lex.Cst}. Two variables are added:
\code{lex.Fail}, an indicator of wheter an event for the transition
names in \code{lex.Tr} has occurred or not. \code{lex.Tr} is a factor
with levels made up of combinations of the levels of \code{lex.Cst} and
\code{lex.Xst} that do occur together in \code{x}, joined by a
"\code{->}".  }
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
ls.dmi <- stack( dmi )
str( ls.dmi )
# Check that all the transitions and person-years got across.
with( ls.dmi, rbind( table(lex.Fail,lex.Tr),
                     tapply(lex.dur,lex.Tr,sum) ) )
}
\seealso{
\code{\link{splitLexis}}
\code{\link{cutLexis}}
\code{\link{Lexis}}
}
\keyword{survival}
