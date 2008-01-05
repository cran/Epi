\name{tab}
\alias{tab}
\alias{tab.Lexis}
\title{
  Tabulate transitions and risk time from a Lexis object
  }
\description{
  A two-way table of records and transitions classified by states
  (\code{lex.Cst} and \code{lex.Xst}), as well the risk time in each state.
  }
\usage{
  \method{tab}{Lexis}( x, simplify=TRUE, scale=1, ... )
  }
\arguments{
  \item{ x }{A Lexis object.}
  \item{ simplify }{Should rows with 0 follow-up time be dropped?}
  \item{scale}{Scaling factor for the rates. The calculated rates are
               multiplied by this number.}
  \item{ ... }{Other parameters - ignored}
}
\details{
  The first part of the table gives the number of rows in the Lexis object,
  classified by states. The \code{D} column is the number of transitions i.e.
  rows where \code{lex.Cst != lex.Xst}, and the final column \code{Y}
  is the person-years spent in each state.
  }
\value{
  A matrix with rows classified by states where persons spend time, and
  columns classified by stated to which persons transit. Two additonal
  columns of transitions from and time spent in each state is added.
  }
\author{Bendix Carstensen, \email{bxc@steno.dk}}
\examples{
data( nickel )
nic <- Lexis( data=nickel,
             entry=list(age=agein),
              exit=list(age=ageout,cal=ageout+dob,tfh=ageout-age1st),
# Lung cancer deaths and other deaths are coded 1 and 2
       exit.status=(icd > 0) + (icd \%in\% c(162,163))  )
str( nic )
head( nic )
tab.Lexis( nic )
tab( nic )
}
\keyword{survival}