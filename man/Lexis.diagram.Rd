\name{Lexis.diagram}
\alias{Lexis.diagram}
\title{Plot a Lexis diagram}
\description{
   Draws a Lexis diagram, optionally with life lines from a cohort, and
   with lifelines of a cohort if supplied. Intended for presentation purposes.
}
\usage{
Lexis.diagram( age = c( 0, 60), 
              alab = "Age",
              date = c( 1940, 2000 ),
              dlab = "Calendar time",
               int = 5,
           lab.int = 2*int, 
          col.life = "black",
          lwd.life = 2,
          age.grid = TRUE,
         date.grid = TRUE,
          coh.grid = FALSE,
          col.grid = gray(0.7),
          lwd.grid = 1,     
               las = 1,
        entry.date = NA,
         entry.age = NA,
         exit.date = NA,
          exit.age = NA,
         risk.time = NA,
        birth.date = NA,
              fail = NA,
          cex.fail = 1.1, 
          pch.fail = c(NA,16),
          col.fail = rep( col.life, 2 ),
              data = NULL, ... )
}
\arguments{
 \item{age}{Numerical vector of length 2, giving the age-range for the diagram}
 \item{alab}{Label on the age-axis.}
 \item{date}{Numerical vector of length 2, giving the calendar
    time-range for the diagram} 
 \item{dlab}{label on the calendar time axis.}
 \item{int}{The interval between grid lines in the diagram. If a
vector of length two is given, the first value will be used for
spacing of age-grid and the second for spacing of the date grid.}
 \item{lab.int}{The interval between labelling of the grids.}
 \item{col.life}{Colour of the life lines.}
 \item{lwd.life}{Width of the life lines.}
 \item{age.grid}{Should grid lines be drawn for age?}
 \item{date.grid}{Should grid lines be drawn for date?}
 \item{coh.grid}{Should grid lines be drawn for birth cohorts (diagonals)?}
 \item{col.grid}{Colour of the grid lines.}
 \item{lwd.grid}{Width of the grid lines.}
 \item{las}{How are the axis labels plotted?}
 \item{entry.date, entry.age, exit.date, exit.age, risk.time,
       birth.date}{Numerical vectors defining lifelines to be plotted
       in the diagram. At least three must be given to produce lines.
       Not all subsets of three will suffice, the given subset has to 
       define life lines.
       If insufficient data is given, no life
       lines are produced.}
 \item{fail}{Logical of event status at exit for the persons whose 
       life lines are plotted.}
 \item{pch.fail}{Symbols at the end of the life lines for censorings
   (\code{fail==0}) and failures (\code{fail != 0}).}
 \item{cex.fail}{Expansion of the status marks at the end of life lines.}
 \item{col.fail}{Character vector of length 2 giving the colour of the
   failure marks for censorings and failures respectively.}
  \item{data}{Dataframe in which to interpret the arguments.}
  \item{...}{Arguments to be passed on to the initial call to plot.}
}
\value{
If sufficient information on lifelines is given, a data frame with
one row per person and columns with entry ages and dates, birth date,
risk time and status filled in.

Side effect: a plot of a Lexis diagram is produced with the life lines
in it is produced. This will be the main reason for using the
function. If the primary aim is to illustrate follow-up of a cohort, then
it is better to represent the follow-up in a \code{\link{Lexis}} object, and
use the generic \code{\link{plot.Lexis}} function.
}
\details{
The default unit for supplied variables are (calendar) years.
If any of the variables \code{entry.date}, \code{exit.date} or
\code{birth.date} are of class "\code{Date}" or if any of the variables
\code{entry.age}, \code{exit.age} or \code{risk.time} are of class
"\code{difftime}", they will be converted to calendar years, and plotted
correctly in the diagram. The returned dataframe will then have colums of
classes "\code{Date}" and  "\code{difftime}".
}
\author{
  Bendix Carstensen,
  \url{http://bendixcarstensen.com}
}
\examples{
Lexis.diagram( entry.age = c(3,30,45),
               risk.time = c(25,5,14),
              birth.date = c(1970,1931,1925.7),
                    fail = c(TRUE,TRUE,FALSE) )
LL <- Lexis.diagram( entry.age = sample( 0:50, 17, replace=TRUE ),
                     risk.time = sample( 5:40, 17, r=TRUE),
                    birth.date = sample( 1910:1980, 17, r=TRUE ),
          		  fail = sample( 0:1, 17, r=TRUE ), 
	              cex.fail = 1.1,
	              lwd.life = 2 )
# Identify the persons' entry and exits
text( LL$exit.date, LL$exit.age, paste(1:nrow(LL)), col="red", font=2, adj=c(0,1) )
text( LL$entry.date, LL$entry.age, paste(1:nrow(LL)), col="blue", font=2, adj=c(1,0) )
data( nickel )
attach( nickel )
LL <- Lexis.diagram( age=c(10,100), date=c(1900,1990), 
	             entry.age=age1st, exit.age=ageout, birth.date=dob, 
		     fail=(icd \%in\% c(162,163)), lwd.life=1,
		     cex.fail=0.8, col.fail=c("green","red") )
abline( v=1934, col="blue" )
nickel[1:10,]
LL[1:10,]
}
\keyword{hplot}
\keyword{dplot}
\seealso{
  \code{\link{Life.lines}},
  \code{\link{Lexis.lines}}
}
