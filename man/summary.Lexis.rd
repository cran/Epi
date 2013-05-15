\name{summary.Lexis}
\alias{summary.Lexis}
\alias{print.summary.Lexis}
\title{
  Summarize transitions and risk time from a Lexis object
  }
\description{
  A two-way table of records and transitions classified by states
  (\code{lex.Cst} and \code{lex.Xst}), as well the risk time in each state.
  }
\usage{
  \method{summary}{Lexis}( object, simplify=TRUE, scale=1, by=NULL, Rates=FALSE, ... )
  \method{print}{summary.Lexis}( x, ..., digits=2 )
  }
\arguments{
  \item{object}{A Lexis object.}
  \item{simplify}{Should rows with 0 follow-up time be dropped?}
  \item{scale}{Scaling factor for the rates. The calculated rates are
    multiplied by this number.}
  \item{by}{Character vector of name(s) of variable(s) in
    \code{object}. Used to give a separate summaries for subsets of
    \code{object}. If longer than than 1, the interaction between that
    variables is used to stratify the summary. It is also possible to
    supply a vector of length \code{nrow(object)}, and the distinct
    values of this will be used to stratify the summary.}
\item{Rates}{Should a component with transition rates be returned (and
  printed) too?}
  \item{x}{A \code{summary.Lexis} object.}
  \item{digits}{How many digits should be used for printing?}
  \item{ ... }{Other parameters - ignored}
}
\value{
  An object of class \code{summary.Lexis}, a list with two components,
  \code{Transitions} and \code{Rates}, each one a matrix with rows
  classified by states where persons spent time, and columns classified
  by states to which persons transit. The \code{Transitions} contains
  number of transitions and has 4 extra columns with number of records,
  total number of events, total risk time and number of person
  contributing attached.  The \code{Rates} contains the transitions
  rates.

  If the argument \code{Rates} is FALSE (the default), then only the
  first component of the list is returned.
}
\author{Bendix Carstensen, \email{bxc@steno.dk}}
\examples{
data( nickel )
# Lung cancer deaths and other deaths are coded 1 and 2
nic <- Lexis( data=nickel,
             entry=list(age=agein),
              exit=list(age=ageout,cal=ageout+dob,tfh=ageout-age1st),
       exit.status=factor( (icd > 0) + (icd \%in\% c(162,163)),
                           labels=c("Alive","Other","Lung") ) )
str( nic )
head( nic )
summary( nic )
# More detailed summary, by exposure level
summary( nic, by=nic$exposure>5, Rates=TRUE, scale=100 )
}
\keyword{survival}