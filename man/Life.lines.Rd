\name{Life.lines}
\alias{Life.lines}
\title{
  Compute dates/ages for life lines in a Lexis diagram
}
\description{
  Fills out the missing information for follow up of persons in a Lexis
  diagram if sufficient information is given.
}
\usage{
Life.lines( entry.date = NA,
             exit.date = NA,
            birth.date = NA,
             entry.age = NA,
              exit.age = NA,
             risk.time = NA )
}
\arguments{
 \item{entry.date, exit.date,birth.date, entry.age, exit.age,
       risk.time}{Vectors defining lifelines to be plotted 
       in the diagram. At least three must be given to produce a result.
       Not all subsets of three will suffice, the given subset has to 
       define life lines. If insufficient data is given, nothing is
       returned and a warning is given.}
}
\value{
  Data frame with variables \code{entry.date}, \code{entry.age},
  \code{exit.date}, \code{exit.age}, \code{risk.time},
  \code{birth.date}, with all entries computed for each person. If any
  of \code{entry.date}, \code{exit.date} or \code{birth.date} are of
  class \code{Date} or if any of \code{entry.age}, \code{exit.age} or
  \code{risk.time} are of class \code{difftime} the date variables will
  be of class  \code{Date} and the other three of class
  \code{difftime}. 
}
\examples{
( Life.lines( entry.age = c(3,30,45),
              risk.time = c(25,5,14),
             birth.date = c(1970,1931,1925.7) ) )

# Draw a Lexis diagram
Lexis.diagram()

# Compute entry and exit age and date.
( LL <-  Life.lines( entry.age = c(3,30,45),
                     risk.time = c(25,5,14),
                    birth.date = c(1970,1931,1925.7) ) )
segments( LL[,1], LL[,2], LL[,3], LL[,4] ) # Plot the life lines.

# Compute entry and exit age and date, supplying a date variable
bd <- ( c(1970,1931,1925.7) - 1970 ) * 365.25
class( bd ) <- "Date"
( Life.lines( entry.age = c(3,30,45),
              risk.time = c(25,5,14),
             birth.date = bd ) )
}
\seealso{
  \code{\link{Lexis.diagram}},
  \code{\link{Lexis.lines}}
}
\keyword{ manip }
\keyword{ dplot }
