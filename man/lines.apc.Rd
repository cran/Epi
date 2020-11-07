\name{apc.lines}
\alias{lines.apc}
\alias{apc.lines}
\title{
  Plot APC-estimates in an APC-frame.
}
\description{
  When an APC-frame has been produced by \code{\link{apc.frame}}, this
  function draws a set of estimates from an APC-fit in the frame. An
  optional drift parameter can be added to the period parameters and
  subtracted from the cohort and age parameters.
}
\usage{
\method{lines}{apc}( x, P, C,
        scale = c("log","ln","rates","inc","RR"),
    frame.par = options()[["apc.frame.par"]],
        drift = 0,
           c0 = median( C[,1] ),
           a0 = median( A[,1] ),
           p0 = c0 + a0,
           ci = rep( FALSE, 3 ),
          lwd = c(3,1,1),
          lty = 1,
          col = "black",
         type = "l",
        knots = FALSE,
        shade = FALSE,
          ... )
 apc.lines( x, P, C,
        scale = c("log","ln","rates","inc","RR"),
    frame.par = options()[["apc.frame.par"]],
        drift = 0,
           c0 = median( C[,1] ),
           a0 = median( A[,1] ),
           p0 = c0 + a0,
           ci = rep( FALSE, 3 ),
          lwd = c(3,1,1),
          lty = 1,
          col = "black",
         type = "l",
        knots = FALSE,
        shade = FALSE,
          ... )
}
\arguments{
  \item{x}{If an \code{apc}-object, (see \code{\link{apc.fit}}), then
    the arguments \code{P}, \code{C}, \code{c0}, \code{a0} and \code{p0}
    are ignored, and the estimates from \code{x} plotted.

    Can also be a 4-column matrix with columns age, age-specific
    rates, lower and upper c.i., in which case period and cohort effects
    are taken from the arguments \code{P} and \code{C}.
    }
  \item{P}{Period effects. Rate-ratios. Same form as for the age-effects.}
  \item{C}{Cohort effects. Rate-ratios. Same form as for the age-effects.}
  \item{scale}{Are effects given on a log-scale? Character variable, one
    of \code{"log"}, \code{"ln"}, \code{"rates"}, \code{"inc"},
    \code{"RR"}. If \code{"log"} or \code{"ln"} it is assumed that
    effects are log(rates) and log(RRs) otherwise the actual effects are
    assumed given in \code{A}, \code{P} and \code{C}. If \code{A} is of
    class \code{apc}, it is assumed to be \code{"rates"}.}
  \item{frame.par}{2-element vector with the cohort-period offset and
    RR multiplicator. This will typically be the result from the call of
    \code{\link{apc.frame}}. See this for details.}
  \item{drift}{The drift parameter to be added to the period effect. If
    \code{scale="log"} this is assumed to be on the log-scale, otherwise
    it is assumed to be a multiplicative factor per unit of the first
    columns of  \code{A}, \code{P} and \code{C} }
  \item{c0}{The cohort where the drift is assumed to be 0; the subtracted
  drift effect is \code{drift*(C[,1]-c0)}.}
  \item{a0}{The age where the drift is assumed to be 0.}
  \item{p0}{The period where the drift is assumed to be 0.}
  \item{ci}{Should confidence interval be drawn. Logical or
    character. If character, any occurrence of \code{"a"} or \code{"A"}
    produces confidence intervals for the age-effect. Similarly for
    period and cohort.}
  \item{lwd}{Line widths for estimates, lower and upper confidence limits.}
  \item{lty}{Linetypes for the three effects.}
  \item{col}{Colours for the three effects.}
  \item{type}{What type of lines / points should be used.}
  \item{knots}{Should knots from the model be shown?}
  \item{shade}{Should confidence intervals be plotted as shaded areas?
    If true, the setting of \code{ci} is ignored.}
  \item{...}{Further parameters to be transmitted to \code{points}
    \code{lines}, \code{matpoints} or \code{matlines} used
    for plotting the three sets of curves.}
}
\details{
  There is no difference between the functions \code{apc.lines} and
    \code{lines.apc}, except the the latter is the \code{lines} method
    for \code{apc} objects.
  
  The drawing of three effects in an APC-frame is a rather trivial task,
  and the main purpose of the utility is to provide a function that
  easily adds the functionality of adding a drift so that several sets
  of lines can be easily produced in the same frame.
}
\value{
  \code{apc.lines} returns (invisibly) a list of three matrices of the
  effects plotted.
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://bendixcarstensen.com}
}
\seealso{
  \code{\link{apc.frame}}, \code{\link{pc.lines}}, \code{\link{apc.fit}}, \code{\link{apc.plot}}
}
\keyword{hplot}
