\name{transform.Lexis}
\alias{transform.Lexis}
\alias{Relevel.Lexis}
\alias{transform.stacked.Lexis}
\alias{factorize}
\alias{factorize.Lexis}
\alias{levels.Lexis}
\alias{order.Lexis}
\alias{sort.Lexis}
\title{Transform a Lexis (or stacked.Lexis) objects}
\description{
  Modify a Lexis object.
}
\usage{
\method{transform}{Lexis}( `_data`, \dots)
\method{Relevel}{Lexis}( x, states, print = TRUE, \dots )
\method{levels}{Lexis}( x )
\method{factorize}{Lexis}( x, states, print = TRUE, \dots )
%\method{order}{Lexis}( x, \dots )
%\method{sort}{Lexis}( x, \dots )
\method{transform}{stacked.Lexis}( `_data`, \dots)
}
\arguments{
  \item{_data}{an object of class \code{Lexis}.}
  \item{x}{an object of class \code{Lexis}.}
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
  method for data frames, but keeps the \code{Lexis} attributes.

  \code{factorize} transforms the variables
  \code{lex.Cst} and \code{lex.Xst} to factors with identical set of
  levels, optionally with names given in \code{states}, and optionally
  collapsing states.

  \code{Relevel} is merely an alias for
  \code{factorize}, since the function does the same as
  \code{\link{Relevel}}, but for both the factors \code{lex.Cst} and
  \code{lex.Xst}. A default sideeffect is to produce a table of old
  states versus new states if \code{states} is a list. Unlike
  \code{Relevel} for factors, \code{Relevel.Lexis} does not accept a
  matrix as a second argument - the number of levels of \code{lex.Cst}
  is rarely (if ever) large.

  If \code{states} is \code{NULL}, as when for example the argument is
  not passed to the function, the returned object have levels of
  \code{lex.Cst}, \code{lex.Xst} (and for \code{stacked.Lexis} objects
  \code{lex.Tr}) shaved down to the actually occurring values.

%  \code{order} returns the order of the rows in a Lexis object to sort
%  it by (\code{lex.id},\code{timeScales[x]}).

%  \code{sort} returns the Lexis object sorted by
%  (\code{lex.id},\code{timeScales[x]}). 
}
\value{
  A transformed \code{Lexis} object.

  The function \code{levels} returns the names of the states (levels of
  the factors \code{lex.Cst} and \code{lex.Xst}.
  }
\author{Martyn Plummer, Bendix Carstensen}
\seealso{\code{\link{Lexis}},
         \code{\link{merge.Lexis}},
         \code{\link{subset.Lexis}},
         \code{\link{subset.stacked.Lexis}},
         \code{\link{Relevel}}}
\examples{
data( nickel )
nic <- Lexis( data = nickel,
                id = id,
             entry = list(age=agein),
              exit = list(age=ageout,cal=ageout+dob,tfh=ageout-age1st),
    ## Lung cancer deaths are coded 2 and other deaths are coded 1
       exit.status = ( (icd > 0) + (icd \%in\% c(162,163)) ) )
str( nic )
levels( nic )
nit <- transform( nic, cumex = exposure*(agein-age1st) )
str( nit )
## It is still a Lexis object!
summary( nic )
nix <- factorize.Lexis( nic, c("Alive","Lung","Dead"))
niw <- factorize.Lexis( nix, c("Alive","Pulm","Mort"))
niz <- factorize.Lexis( niw, states=list("Alive",c("Pulm","Mort")), coll=" \n& ")
boxes( niw, boxpos=TRUE )
par( new=TRUE )
boxes( niz, boxpos=TRUE )
siw <- stack( niw )
str( siw )
}
\keyword{manip}
