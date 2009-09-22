\name{transform.Lexis}
\alias{transform.Lexis}
\title{Transform a Lexis objects}
\description{
  Transform a Lexis object
}
\usage{
%transform(`_data`, \dots)
\method{transform}{Lexis}(`_data`, \dots)
}
\arguments{
  \item{_data}{an object of class \code{Lexis}.}
  \item{\dots}{additional arguments to be passed to \code{transform.data.frame}.}
}
\details{
  The transform method for \code{Lexis} objects works exactly as the
  method for data frames.
}
\value{
  A transformed \code{Lexis} object.
}
\author{Martyn Plummer}
\seealso{\code{\link{Lexis}},
         \code{\link{merge.Lexis}},
         \code{\link{subset.Lexis}}}
\keyword{manip}
