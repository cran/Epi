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
  \method{summary}{Lexis}( object, simplify=TRUE, scale=1, ... )
  \method{print}{summary.Lexis}( x, ..., digits=2 )
  }
\arguments{
  \item{object}{A Lexis object.}
  \item{x}{A \code{summary.Lexis} object.}
  \item{simplify}{Should rows with 0 follow-up time be dropped?}
  \item{scale}{Scaling factor for the rates. The calculated rates are
               multiplied by this number.}
  \item{digits}{How many digits should be used for printing?}
  \item{ ... }{Other parameters - ignored}
}
\value{
  An object of class \code{summary.Lexis}, a list with two components,
  \code{Transitions} and \code{Rates}, each one
  a matrix with rows classified by states where persons spend time, and
  columns classified by stated to which persons transit. The \code{Transitions}
  contains number of transitions and
  has two extra columns of total number events and total risk time attached.
  The \code{Rates} contians the transitions rates.
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
}
\keyword{survival}