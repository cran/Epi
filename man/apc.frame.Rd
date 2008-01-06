\name{apc.frame}
\alias{apc.frame}
\title{
  Produce an empty frame for display of parameter-estimates from
  Age-Period-Cohort-models.
}
\description{
  
  A plot is generated where both the age-scale and the cohort/period
  scale is on the x-axis. The left vertical axis will be a logarithmic
  rate scale referring to age-effects and the right a logarithmic
  rate-ratio scale of the same relative extent as the left referring to
  the cohort and period effects (rate ratios).

  Only an empty plot frame is generated. Curves or points must be added
  with \code{points}, \code{lines} or the special utility function
  \code{\link{apc.lines}}. 
}
\usage{
  apc.frame( a.lab,
            cp.lab,
             r.lab,
            rr.lab = r.lab / rr.ref,
            rr.ref = r.lab[length(r.lab)/2],
             a.tic = a.lab,
            cp.tic = cp.lab,
             r.tic = r.lab,
            rr.tic = r.tic / rr.ref,
           tic.fac = 1.3,
             a.txt = "Age",
            cp.txt = "Calendar time",
             r.txt = "Rate per 100,000 person-years",
            rr.txt = "Rate ratio",
               gap = diff(range(c(a.lab, a.tic)))/3,
          col.grid = gray(0.85),
             sides = c(1,2,4) )
}
\arguments{ 
  \item{a.lab}{Numerical vector of labels for the age-axis.}
  \item{cp.lab}{Numerical vector of labels for the cohort-period axis.}
  \item{r.lab}{Numerical vector of labels for the rate-axis (left vertical)}
  \item{rr.lab}{Numerical vector of labels for the RR-axis (right vertical)}
  \item{rr.ref}{At what level of the rate scale is the RR=1 to be.}
  \item{a.tic}{Location of additional tick marks on the age-scale}
  \item{cp.tic}{Location of additional tick marks on the cohort-period-scale}
  \item{r.tic}{Location of additional tick marks on the rate-scale}
  \item{rr.tic}{Location of additional tick marks on the RR-axis.}
  \item{tic.fac}{Factor with which to diminish intermediate tick marks}
  \item{a.txt}{Text for the age-axis (left part of horizontal axis).}
  \item{cp.txt}{Text for the cohort/period axis (right part of
    horizontal axis).}
  \item{r.txt}{Text for the rate axis (left vertical axis).}
  \item{rr.txt}{Text for the rate-ratio axis (right vertical axis)}
  \item{gap}{Gap between the age-scale and the cohort-period scale}
  \item{col.grid}{Colour of the grid put in the plot.}
  \item{sides}{Numerical vector indicating on which sides axes should
    be drawn and annotated. This option is aimed for multi-panel
    displays where axes only are put on the outer plots.} 
}
\details{
  The function produces an empty plot frame for display of results
  from an age-period-cohort model, with age-specific rates in the left
  side of the frame and cohort and period rate-ratio parameters in the
  right side of the frame. There is a gap of \code{gap} between the
  age-axis and the calendar time axis, vertical grid lines at
  \code{c(a.lab,a.tic,cp.lab,cp.tic)}, and horizontal grid lines at
  \code{c(r.lab,r.tic)}.
  
  The function returns a numerical vector of
  length 2, with names \code{c("cp.offset","RR.fac")}. The y-axis for
  the plot will be a rate scale for the age-effects, and the x-axis will
  be the age-scale. The cohort and period effects are plotted by
  subtracting the first element (named \code{"cp.offset"}) of the returned result
  form the cohort/period, and multiplying  the rate-ratios by the second
  element of the returned result (named \code{"RR.fac"}).
}
\value{
  A numerical vector of length two, with names
  \code{c("cp.offset","RR.fac")}. The first is the offset for the cohort
  period-axis, the second the multiplication factor for the rate-ratio
  scale.

  Side-effect: A plot with axes and grid lines but no points or curves.
}
\references{
  \url{http://www.pubhealth.ku.dk/~bxc/APC/notes.pdf}
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://www.pubhealth.ku.dk/~bxc/}
} 
\examples{
par( mar=c(4,4,1,4) )
res <-
apc.frame( a.lab=seq(30,90,20), cp.lab=seq(1880,2000,30), r.lab=c(1,2,5,10,20,50),
           a.tic=seq(30,90,10), cp.tic=seq(1880,2000,10), r.tic=c(1:10,1:5*10),
           gap=27 )
res
# What are the axes actually?
par(c("usr","xlog","ylog"))
# How to plot in the age-part: a point at (50,10)
points( 50, 10, pch=16, cex=2, col="blue" )
# How to plot in the cohort-period-part: a point at (1960,0.3)
points( 1960-res[1], 0.3*res[2], pch=16, cex=2, col="red" )
}
\seealso{
  \code{\link{apc.lines},\link{apc.fit}}
}
\keyword{hplot}
