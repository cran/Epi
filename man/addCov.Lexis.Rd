\name{addCov.Lexis}
\alias{addCov.Lexis}
\title{
  Add covariates (typically clinical measurements) taken at known times
  to a Lexis object. 
}
\description{
  When follow-up in a multistate model is represented in a
  \code{\link{Lexis}} object we may want to add information on
  covariates, for example clinical measurements, obtained at different
  times. This function cuts the follow-up time (see
  \code{\link{cutLexis}}) at the times of measurement and carries the
  measurements forward in time to the next measurement occasion.
}
\usage{
\method{addCov}{Lexis}( Lx,
                      clin,
                 timescale = 1,
                     exnam,
                       tfc = "tfc",
                 addScales = FALSE )
}
\arguments{
  \item{Lx}{
    A Lexis object with follow-up of a cohort.  
    }
  \item{clin}{
    A data frame with covariates to add (typically clinical
    measurements). Must contain a variable \code{lex.id} identifying the
    persons represented in \code{Lx}, as well as a variable with the
    same name as one of the \code{\link{timeScales}} in \code{Lx},
    identifying the time at which covariates are measured.

    The times must be unique within each person; if not records with
    duplicate times are discarded, and a warning issued. This is done
    using \code{duplicated}, so not very well-defined, it is advisable
    that you do this by yourself.
    }
  \item{timescale}{
    Numerical or character. Number or name of a timescale in
    \code{Lx}. The \code{clin} data frame must have a variable of this
    name indicating the time at which the covariate measurements were
    taken.
    }
    \item{exnam}{
    Character. Name of the variable in \code{clin} with the examination
    names (such as \code{wave1}, \code{wave2} etc.). Values may not be
    repeated within person. Will be carried over to the resulting
    \code{Lexis} object. If there is no such variable in \code{clin} it
    will be constructed; if the argument is omitted, a variable called
    \code{exnam} with values \code{ex1}, \code{ex2} etc. will be
    constructed.
    }
  \item{tfc}{      
    Character (\code{t}ime \code{f}rom \code{c}ovariate). Name of the
    variable in the result which will contain the time since the most
    recent covariate date. If the argument is omitted a variable called
    \code{tfc} will be constructed. If \code{addScales} is \code{TRUE}
    this is included among the time scales.
    }
  \item{addScales}{
    Logical. Should timescales representing time since each examination time
    time be added? They will be named \code{paste("tf",exnam)}.
    }
}
\value{
    A \code{Lexis} object representing the same follow-up as \code{Lx}, 
    with cuts added at the times of examination, and covariate
    measurements added for all records representing follow-up after the
    most recent time of measurement. 
}
\details{
    The implementation is clumpy, the function is slow.
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

# Data frame with clinical examination data, date of examination in per
clin <- data.frame( lex.id = c(1,1,3,2),
                       per = c(1977.3,1971.7,1996.2,1990.6),
                        bp = c(120,140,160,157),
                      chol = c(5,7,8,9),
                      xnam = c("X2","X1","X1","X2") )
Lx
clin 

# Different behavours using exnam and addScales
addCov.Lexis( Lx, clin )
addCov.Lexis( Lx, clin, exnam="xnam" )
addCov.Lexis( Lx, clin, addScales=TRUE )
addCov.Lexis( Lx, clin, addScales=TRUE, exnam="xnam" )

# Works with time split BEFORE
Lb <- addCov.Lexis( splitLexis( Lx,
                                time.scale="age",
                                breaks=seq(0,80,5) ),
                    clin,
                    exnam="clX" )
Lb
# ...or AFTER, but only if addScales is given
La <- splitLexis( addCov.Lexis( Lx,
                              clin,
                             exnam = "xnam" ),
                  breaks=seq(0,80,5),
                  time.scale="age" )
La
La <- splitLexis( addCov.Lexis( Lx,
                              clin,
                             exnam = "xnam",
                         addScales = TRUE ),
                  breaks=seq(0,80,5),
                  time.scale="age" )
La
}
\keyword{ survival }
