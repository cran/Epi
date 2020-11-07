\name{Lexis.lines}
\alias{Lexis.lines}
\title{Draw life lines in a Lexis diagram.}
\description{
  Add life lines to a Lexis diagram.
}
\usage{
Lexis.lines( entry.date = NA,
              exit.date = NA,
             birth.date = NA,
              entry.age = NA,
               exit.age = NA,
              risk.time = NA,
               col.life = "black",
               lwd.life = 2,
                   fail = NA,
               cex.fail = 1,
               pch.fail = c(NA, 16),
               col.fail = col.life,
                   data = NULL )
}
\arguments{
 \item{entry.date, entry.age, exit.date, exit.age, risk.time,
       birth.date}{Numerical vectors defining lifelines to be plotted
       in the diagram. At least three must be given to produce lines.
       Not all subsets of three will suffice, the given subset has to 
       define life lines. If insufficient data is given, no life
       lines are produced.}
 \item{col.life}{Colour of the life lines.}
 \item{lwd.life}{Width of the life lines.}
 \item{fail}{Logical of event status at exit for the persons whose 
       life lines are plotted.}
 \item{cex.fail}{The size of the status marks at the end of life lines.}
 \item{pch.fail}{The status marks at the end of the life lines.}
 \item{col.fail}{Colour of the marks for censorings and failures
   respectively.} 
 \item{data}{Data frame in which to interpret values.}
}
\value{
  If sufficient information on lifelines is given, a data frame with one
  row per person and columns with entry ages and dates, birth date, risk
  time and status filled in.

  Side effect: Life lines are added to an existing Lexis
  diagram. Lexis.lines adds life lines to an existing plot. 
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://bendixcarstensen.com}
}
\examples{
Lexis.diagram( entry.age = c(3,30,45),
               risk.time = c(25,5,14),
              birth.date = c(1970,1931,1925.7),
                    fail = c(TRUE,TRUE,FALSE) )
Lexis.lines( entry.age = sample( 0:50, 100, replace=TRUE ),
             risk.time = sample( 5:40, 100, r=TRUE),
            birth.date = sample( 1910:1980, 100, r=TRUE ),
                  fail = sample(0:1,100,r=TRUE),
              cex.fail = 0.5,
              lwd.life = 1 )
}
\keyword{ hplot }
\keyword{ dplot }
\seealso{
  \code{\link{Lexis.diagram}},
  \code{\link{Life.lines}}
}
