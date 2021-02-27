\name{transform.Lexis}
\alias{transform.Lexis}
\alias{Relevel.Lexis}
\alias{transform.stacked.Lexis}
\alias{factorize}
\alias{factorize.Lexis}
\alias{levels.Lexis}
\alias{order.Lexis}
\alias{orderLexis}
\alias{sortLexis}
\title{Transform a Lexis (or stacked.Lexis) object}
\description{
  Modify a Lexis object.
}
\usage{
\method{factorize}{Lexis}(x, ...)
\method{Relevel}{Lexis}(x, ref, \dots)
\method{levels}{Lexis}(x)
\method{transform}{Lexis}(`_data`, \dots)
\method{transform}{stacked.Lexis}(`_data`, \dots)
 order.Lexis(x)
  orderLexis(x)
   sortLexis(x)
}
\arguments{
  \item{_data}{an object of class \code{Lexis}.}
  \item{x}{an object of class \code{Lexis}.}
  \item{ref}{New names (or order) of the factor levels (states) for
    \code{lex.Cst} 
    and \code{lex.Xst}. Can be a list, in which case some levels are
    collapsed, see the documentation for \code{\link{Relevel}}. No
    sanity check for the latter type of operation is undertaken.}
  \item{\dots}{Additional arguments to be passed to
    \code{\link{transform.data.frame}}, \code{\link{Relevel.factor}}.}
%   \code{order} or \code{sort}.}
}
\details{
  The transform method for \code{Lexis} objects works exactly as the
  method for data frames, but keeps the \code{Lexis} attributes.

  \code{factorize} transforms the variables
  \code{lex.Cst} and \code{lex.Xst} to factors with identical sets of
  levels.
  
  \code{Relevel} does the same as \code{\link{Relevel.factor}}, but for
  both the factors \code{lex.Cst} and \code{lex.Xst} in \code{x}.
  \code{lex.Cst} and \code{lex.Xst} must be factors with the same
  levels. They can be made so by \code{factorize}.

  If \code{ref} is an integer or character vector, the levels of
  \code{lex.Cst} and \code{lex.Xst} are permuted to match the order of
  \code{ref}. 
  
  If \code{ref} is \code{NULL}, as when for example the argument is
  not passed to the function, the returned object have levels of
  \code{lex.Cst}, \code{lex.Xst} (and for \code{stacked.Lexis} objects
  \code{lex.Tr}) shaved down to the actually occurring values; that is,
  empty levels are discarded.

  \code{order.Lexis} returns the order of the rows in a Lexis object to sort
  it by ()\code{lex.id},\code{ts}), where \code{ts} is a timescale in
  the Lexis object with no \code{NA}s. \code{orderLexis} is just a synonym.

  \code{sortLexis} returns the Lexis object sorted by
  (\code{lex.id}, \code{ts}) where \code{ts} is one of the
  \code{\link{timeScales}} with no \code{NA}s.   
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
         \code{\link{Relevel}},
         \code{\link{transient}},
         \code{\link{absorbing}}}
\examples{
data( nickel )
nic <- Lexis( data = nickel,
                id = id,
             entry = list(age = agein),
              exit = list(age = ageout,
                          cal = ageout+dob,
                          tfh = ageout-age1st),
# Lung cancer deaths end as 2 and other deaths as 1
       exit.status = factor((icd > 0) + (icd \%in\% c(162,163)),
                            labels = c("Alive","Dead","Lung") ) )
str( nic )
levels( nic )
nit <- transform( nic, cumex = exposure * (agein - age1st) )
str( nit )

# It is still a Lexis object!
summary(nic)

# change order of levels
nix <- Relevel(nic, c("Alive", "Lung", "Dead"))
summary(nix)

# change names of levels 
niw <- Relevel(nix, list("Alive" = 1, "Pulm" = "Lung", "Mort" = "Dead"))
summary(niw)
boxes(niw, boxpos = TRUE)

# combine levels
niz <- Relevel(niw, list("Alive", c("Pulm", "Mort")), coll=" \n& ")
summary(niz)
par( new = TRUE )
boxes(niz, boxpos = TRUE)

#stack Lexis object
siw <- stack(niw)
str(siw)
}
\keyword{manip}
