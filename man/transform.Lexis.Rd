\name{transform.Lexis}
\alias{transform.Lexis}
\alias{factorize}
\alias{factorize.Lexis}
\title{Transform a Lexis objects}
\description{
  Modify a Lexis object.
}
\usage{
%transform(`_data`, \dots)
\method{transform}{Lexis}( `_data`, \dots)
\method{factorize}{Lexis}( obj, states, print = TRUE, \dots )
}
\arguments{
  \item{_data}{an object of class \code{Lexis}.}
  \item{obj}{an object of class \code{Lexis}.}
  \item{states}{Names of the factor levels (states) for \code{lex.Cst}
    and \code{lex.Xst}. Can be a list, in which case some levels are
    collapsed, see the documentation for \code{\link{Relevel}}. No
    sanity check for the latter operation is undertaken.}
  \item{print}{Should a conversion between old and new levels be printed?}
  \item{\dots}{Additional arguments to be passed to
    \code{\link{transform.data.frame}} or \code{\link{Relevel}}.}
}
\details{
  The transform method for \code{Lexis} objects works exactly as the
  method for data frames. \code{factorize} transforms the variables
  \code{lex.Cst} and \code{lex.Xst} to factors with identical set of
  levels, optionally with names given in \code{states}, and optionally
  collapsing states. The default is to produce a table of old states
  versus new states if \code{states} is a list.
  }
\value{
  A transformed \code{Lexis} object.
  }
\author{Martyn Plummer, Bendix Carstensen}
\seealso{\code{\link{Lexis}},
         \code{\link{merge.Lexis}},
         \code{\link{subset.Lexis}}}
\examples{
data( nickel )
nic <- Lexis( data = nickel,
                id = id,
             entry = list(age=agein),
              exit = list(age=ageout,cal=ageout+dob,tfh=ageout-age1st),
    ## Lung cancer deaths are coded 2 and other deaths are coded 1
       exit.status = ( (icd > 0) + (icd \%in\% c(162,163)) ) )
str( nic )
nit <- transform( nic, cumex = exposure*(agein-age1st) )
str( nit )
## It is still a Lexis object!
summary( nic )
nix <- factorize.Lexis( nic, c("Alive","Lung","Dead"))
niw <- factorize.Lexis( nix, c("Alive","Pulm","Mort"))
niz <- factorize.Lexis( nic, states=list("Alive",c("Pulm","Mort")), coll=" \n& ")
boxes( niw, boxpos=TRUE )
par( new=TRUE )
boxes( niz, boxpos=TRUE )
}
\keyword{manip}
