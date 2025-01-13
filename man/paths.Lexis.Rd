\name{paths.Lexis}
\alias{paths.Lexis}
\title{
Generate paths travelled through a Lexis multistate model data frame. 
}
\description{
Paths visited in a \code{Lexis} multistate model.
}
\usage{
\method{paths}{Lexis}(Lx, dfr = FALSE, \dots)
}
\arguments{
  \item{Lx}{A \code{\link{Lexis}} object}
  \item{dfr}{Logical. Should results be returned as a data frame with
    columns \code{lex.id} and \code{path}?}
  \item{\dots}{Arguments passed on. Ignored}
}
\value{A factor with levels describing each person's path through states. It is
  of length \code{length(\link{nid}(Lx))}, named by the (character)
  values of \code{Lx$lex.id}. If \code{dfr} is \code{TRUE} a two-column
  data frame is returned.
}
\author{
  Bendix Carstensen, \email{b@bxc.dk}, \url{http://bendixcarstensen.com}
}
\seealso{
\code{\link{cutLexis}},
\code{\link{mcutLexis}},
\code{\link{rcutLexis}},
\code{\link{nid}},
\code{\link{Lexis}}
}
\examples{
# a simple example
example(DMlate)
summary(dmi)
str(paths.Lexis(dmi, dfr = TRUE))
str(pathD <- paths.Lexis(dmi))
cbind(addmargins(table(pathD)))
#
# an example with recurring events
example(steno2)
summary(L4)
str(pathS <- paths.Lexis(L4))
cbind(addmargins(table(pathS)))
  }
\concept{Lexis}
\concept{Multistate model}
\concept{Multistate path}
