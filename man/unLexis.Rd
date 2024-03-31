\name{unLexis}
\alias{unLexis}
\title{Remove Lexis attributes from a \code{\link{Lexis}} object.}
\description{Removes the Lexis attributes, including the class
  \code{Lexis} from a Lexis object.}
\usage{unLexis(Lx)}
\arguments{\item{Lx}{A Lexis object}}
\value{The input object with "\code{Lexis}" removed from the class attribute.}
\author{Bendix Carstensen}
\seealso{\code{\link{Lexis}}}
\examples{
# A small bogus cohort
xcoh <- structure(list( id = c("A", "B", "C"),
                     birth = c("14/07/1952", "01/04/1954", "10/06/1987"),
                     entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                      exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                      fail = c(1, 0, 1) ),
                    .Names = c("id", "birth", "entry", "exit", "fail"),
                 row.names = c("1", "2", "3"),
                     class = "data.frame")

# Convert the character dates into numerical variables (fractional years)
  xcoh <- cal.yr(xcoh, format="\%d/\%m/\%Y", wh=2:4)
# xcoh <- cal.yr(xcoh, format="%d/%m/%Y", wh=2:4)

# Define a Lexis object with timescales calendar time and age
Lcoh <- Lexis(entry = list(per = entry ),
               exit = list(per = exit,
                           age = exit - birth),
        exit.status = fail,
               data = xcoh)
summary(Lcoh)
try(summary(unLexis(Lcoh)))
}
\concept{data manipulation}

