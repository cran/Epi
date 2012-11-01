\name{stattable.funs}
\alias{count}
\alias{percent}
\alias{ratio}
\title{Special functions for use in stat.table}
\description{
  These functions may be used as \code{contents} arguments to the
  function \code{stat.table}. They are defined internally in
  \code{stat.table} and have no independent existence.
}
\usage{
count(id)
ratio(d,y,scale=1, na.rm=TRUE)
percent(...)
}
\arguments{
  \item{id}{numeric vector in which identical values identify the
    same individual.}
  \item{d, y}{numeric vectors of equal length (\code{d} for Deaths,
    \code{y} for person-Years)}
  \item{scale}{a scalar giving a value by which the ratio should be
    multiplied}
  \item{na.rm}{a logical value indicating whether \code{NA} values should
    be stripped before computation proceeds.}
  \item{...}{a list of variables taken from the \code{index} argument
    to \code{\link{stat.table}}}
}
\value{
  When used as a \code{contents} argument to \code{stat.table}, these
  functions create the following tables:
  \item{\code{count}}{If given without argument (\code{count()}) it
    returns a contingency table of counts. If given an \code{id}
    argument it returns a table of the number of different values of
    \code{id} in each cell, i.e. how many persons contribute in each
    cell.}
  \item{\code{ratio}}{returns a table of values
    \code{scale * sum(d)/sum(y)}}
  \item{\code{percent}}{returns a table of percentages of the
    classifying variables. Variables that are in the \code{index}
    argument to \code{stat.table} but not in the call to
    \code{percent} are used to define strata, within which the
    percentages add up to 100.}
}
\author{Martyn Plummer}
\seealso{\code{\link{stat.table}}}
\keyword{iteration}
\keyword{category}
