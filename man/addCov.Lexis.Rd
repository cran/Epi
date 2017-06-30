\name{addCov.Lexis}
\alias{addCov.Lexis}
\title{
  Add covariates (typically clinical measurements) taken at known times
  to a Lexis object. 
}
\description{
  When follow-up in a multistate model is represented in a
  \code{\link{Lexis}} object we may want to add information on
  covariates, for example clinical measurements, taken at different
  times. This function cuts the follow-up time (see
  \code{\link{cutLexis}}) at the times of measurement and carries the
  measurements forward in time. Several measurements per person is
  allowed in which case measurements from one occasion are only carried
  forward to the next.
}
\usage{
\method{addCov}{Lexis}(Lx, clin,
                       timeScale = 1,
                       addScales = FALSE,
                       exnam,
                       tfc = "tfc")
}
\arguments{
  \item{Lx}{
    A Lexis object with follow-up of a cohort  
    }
  \item{clin}{
    A data frame with the covariates to add (typically clinical
    measurements). Must contain a variable \code{lex.id} identifying the
    persons represented in \code{Lx}, as well as a variable with the
    same name as one of the \code{\link{timeScales}} in \code{Lx},
    identifying the time at which covariates are measured.
    }
  \item{timeScale}{
    Numerical or character. If numerical, the number of a timescale in
    \code{Lx}. The \code{clin} data frame must have a variable of this
    name indicating the time at which the clinical measurements are taken
    (assumed to be all other variables than \code{lex.id},
    \code{timeScale} and \code{exnam}.)
    }
  \item{addScales}{
    Logical. Should timescales representing time since each covariate
    time be added? They will be named \code{paste("tf",exnam)}. Not
    implemented, argument ignored.  
    }
  \item{exnam}{      
    Character. Name of the variable in \code{clin} with the examination
    name (such as \code{wave1}, \code{wave2} etc.). Must not be repeated
    within person. Will be carried over to the resulting \code{Lexis}
    object. If omitted, a variable called \code{exnam} with values
    \code{ex.1}, \code{ex.2} etc. will be constructed.
    }
  \item{tfc}{      
    Character (\code{t}ime \code{f}rom \code{c}ovariate). Name of the
    variable in the result which will contain the time since the most
    recent covariate date. This is not a time scale as it is reset to
    0 at each new covariate time. Also note that this is meaningless if
    you \code{splitLexis} \emph{after} using \code{addCov.Lexis}.
    }
}
\value{
  A \code{Lexis} object representing the same follow-up as in \code{Lx},
  with cuts added at the times of examination, and clinical measurements added
  for all records representing follow-up after the most recent time of
  measurement. 
}
\author{
  Bendix Carstensen, \email{b@bxc.dk}, \url{http://BendixCarstensen.com}
}
\seealso{
\code{\link{cutLexis}},
\code{\link{mcutLexis}},
\code{\link{splitLexis}},
\code{\link{Lexis}}
}
\examples{
# A small bogus cohort
xcoh <- structure( list( id = c("A", "B", "C"),
                      birth = c("1952-07-14", "1954-04-01", "1987-06-10"),
                      entry = c("1965-08-04", "1972-09-08", "1991-12-23"),
                       exit = c("1997-06-27", "1995-05-23", "1998-07-24"),
                       fail = c(1, 0, 1) ),
                     .Names = c("id", "birth", "entry", "exit", "fail"),
                  row.names = c("1", "2", "3"),
                      class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh$bt <- cal.yr( xcoh$birth )
xcoh$en <- cal.yr( xcoh$entry )
xcoh$ex <- cal.yr( xcoh$exit  )

# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( per=en ),
                exit = list( per=ex, age=ex-bt ),
         exit.status = factor( fail, 0:1, c("Alive","Dead") ),
                data = xcoh )
str( Lcoh )
Lx <- Lcoh[,1:7]

# Data frame with clinical examination data, date of examination in
clin <- data.frame( lex.id = c(1,1,3,2),
                       per = c(1977,1971,1996,1990),
                        bp = c(120,140,160,157),
                      chol = c(5,7,8,9) )
Lx
clin 

# Works with time split BEFORE adding clinical data:
Lb <- addCov.Lexis( splitLexis( Lx,
                                time.scale="age",
                                breaks=seq(0,80,5) ),
                    clin,
                    exnam="clX" )
Lb
# With time split AFTER adding clinincal data, variable tfc is meaningless:
La <- splitLexis( addCov.Lexis( Lx,
                                clin,
                                exnam="clX" ),
                  breaks=seq(0,80,5),
                  time.scale="age" )
La
}
\keyword{ survival }
