\name{list.lex}
\alias{list.lex}
\title{Compact printing of Lexis object
}
\description{Often the number of digits after the decimal point ant not
  the sugnificant digits is relevant in printing timescales.
}
\usage{
list.lex( Lx, dig, sh="_" )
}
\arguments{
  \item{Lx}{A \code{\link{Lexis}} object.}
  \item{dig}{Scalar. How many digits should be printed after the decimal
    separator for the numerical variables.}
  \item{sh}{Replacement string for \code{lex.} in names when printing.}
}
\details{The function is designed to make a compact listing of a Lexis
  obejct. Variables are ordered with \code{lex.id}, time scales,
  \code{lex.dur} and state variables first. In order to avoid save
  between columns the \code{lex.} prefix is replaced by an underscore to
  shorten the names.
}
\value{NULL, the functions is used for its compact printing of a Lexis object. 
}
\author{Bendix Carstensen
}
\examples{
# A small bogus cohort
xcoh <- structure( list( id = c("A", "B", "C"),
                      birth = c("14/07/1952", "01/04/1954", "10/06/1987"),
                      entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                       exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                       fail = c(1, 0, 1) ),
                     .Names = c("id", "birth", "entry", "exit", "fail"),
                  row.names = c("1", "2", "3"),
                      class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh <- cal.yr( xcoh, format="\%d/\%m/\%Y", wh=2:4 )
# Make a Lexis obejct
Lcoh <- Lexis( entry = list( per=entry ),
                exit = list( per=exit,
                             age=exit-birth ),
         exit.status = fail,
                data = xcoh )
Lcoh
list.lex( Lcoh, 2 )
}
\keyword{manip}
