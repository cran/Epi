\name{pc.lines}
\alias{pc.points}
\alias{pc.lines}
\alias{pc.matpoints}
\alias{pc.matlines}
\alias{pc.matshade}
\alias{cp.points}
\alias{cp.lines}
\alias{cp.matpoints}
\alias{cp.matlines}
\alias{cp.matshade}
\title{
  Plot period and cohort effects in an APC-frame.
}
\description{
  When an APC-frame has been produced by \code{\link{apc.frame}}, this
  function draws curves or points in the period/cohort part of the frame.
}
\usage{
  pc.points( x, y, ... )
  pc.lines( x, y, ... )
  pc.matpoints( x, y, ... )
  pc.matlines( x, y, ... )
  pc.matshade( x, y, ... )
  cp.points( x, y, ... )
  cp.lines( x, y, ... )
  cp.matpoints( x, y, ... )
  cp.matlines( x, y, ... )
  cp.matshade( x, y, ... )
}
\arguments{
  \item{x}{vector of \code{x}-coordinates.}
  \item{y}{vector or matrix of \code{y}-coordinates.}
  \item{...}{Further parameters to be transmitted to points, lines,
    matpoints, matlines or matshade used for plotting the three sets of curves.}
}
\details{Since the Age-part of the frame is referred to by its real
  coordinates plotting in the calendar time part requires translation
  and scaling to put things correctly there, that is done by the
  functions \code{pc.points} etc.

  The functions \code{cp.points} etc. are just synonyms for these, in
  recognition of the fact that you can never remember whether it is "pc"
  or "cp".
  }
\value{
  The functions return nothing.
  }
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://BendixCarstensen.com}
  }
\seealso{
  \code{\link{apc.frame}}, \code{\link{apc.fit}}, \code{\link{plot.apc}}, \code{\link{lines.apc}}
}
\keyword{hplot}
