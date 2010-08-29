\name{transform.Lexis}
\alias{transform.Lexis}
\alias{factorize}
\alias{factorize.Lexis}
\title{Transform a Lexis objects}
\description{
  Transform a Lexis object
}
\usage{
%transform(`_data`, \dots)
\method{transform}{Lexis}(`_data`, \dots)
\method{factorize}{Lexis}(obj, \dots)
}
\arguments{
  \item{_data}{an object of class \code{Lexis}.}
  \item{obj}{an object of class \code{Lexis}.}
  \item{\dots}{additional arguments to be passed to \code{transform.data.frame}.}
}
\details{
  The transform method for \code{Lexis} objects works exactly as the
  method for data frames.
  
  \code{factorize} simply transforms the variables \code{lex.Cst} and
  \code{lex.Xst} to factors with the same set of levels.
}
\value{
  A transformed \code{Lexis} object.
}
\author{Martyn Plummer, Bendix Carstensen}
\seealso{\code{\link{Lexis}},
         \code{\link{merge.Lexis}},
         \code{\link{subset.Lexis}}}
\keyword{manip}
