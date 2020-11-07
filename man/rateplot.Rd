\name{rateplot}
\alias{rateplot}
\alias{Aplot}
\alias{Pplot}
\alias{Cplot}
\title{
  Functions to plot rates from a table classified by age and
  calendar time (period)
}
\description{  
  Produces plots of rates versus age, connected within period or cohort
  (\code{Aplot}), rates versus period connected within age-groups
  (\code{Pplot}) and rates and rates versus date of birth cohort
  (\code{Cplot}). \code{rateplot} is a wrapper for these, allowing
  to produce the four classical displays with a single call. 
}
\usage{
rateplot( rates,
          which = c("ap","ac","pa","ca"),
            age = as.numeric( dimnames( rates )[[1]] ),
            per = as.numeric( dimnames( rates )[[2]] ),
           grid = FALSE,
         a.grid = grid,
         p.grid = grid,
         c.grid = grid,
          ygrid = grid,
       col.grid = gray( 0.9 ),
          a.lim = range( age, na.rm=TRUE ) + c(0,diff( range( age ) )/30),
          p.lim = range( per, na.rm=TRUE ) + c(0,diff( range( age ) )/30),
          c.lim = NULL, 
           ylim = range( rates[rates>0], na.rm=TRUE ),
             at = NULL,
         labels = paste( at ),
          a.lab = "Age at diagnosis",
          p.lab = "Date of diagnosis",
          c.lab = "Date of birth",
           ylab = "Rates",
           type = "l",
            lwd = 2,
            lty = 1,
         log.ax = "y",
            las = 1,
            ann = FALSE,
          a.ann = ann,
          p.ann = ann,
          c.ann = ann,
          xannx = 1/20,
        cex.ann = 0.8,
         a.thin = seq( 1, length( age ), 2 ),
         p.thin = seq( 1, length( per ), 2 ),
         c.thin = seq( 2, length( age ) + length( per ) - 1, 2 ),
            col = par( "fg" ),
          a.col = col,
          p.col = col,
          c.col = col,
            ... )

Aplot( rates, age = as.numeric( dimnames( rates )[[1]] ),
       per = as.numeric( dimnames( rates )[[2]] ), grid = FALSE,
       a.grid = grid, ygrid = grid, col.grid = gray( 0.9 ),
       a.lim = range( age, na.rm=TRUE ), ylim = range( rates[rates>0], na.rm=TRUE ),
       at = NULL, labels = paste( at ), a.lab = names( dimnames( rates ) )[1],
       ylab = deparse( substitute( rates ) ), type = "l", lwd = 2, lty = 1,
       col = par( "fg" ), log.ax = "y", las = 1, c.col = col, p.col = col,
       c.ann = FALSE, p.ann = FALSE, xannx = 1/20, cex.ann = 0.8,
       c.thin = seq( 2, length( age ) + length( per ) - 1, 2 ),
       p.thin = seq( 1, length( per ), 2 ), p.lines = TRUE,
       c.lines = !p.lines, ... )

Pplot( rates, age = as.numeric( dimnames( rates )[[1]] ),
       per = as.numeric( dimnames( rates )[[2]] ), grid = FALSE,
       p.grid = grid, ygrid = grid, col.grid = gray( 0.9 ),
       p.lim = range( per, na.rm=TRUE ) + c(0,diff(range(per))/30),
       ylim = range( rates[rates>0], na.rm=TRUE ), p.lab = names( dimnames( rates ) )[2],
       ylab = deparse( substitute( rates ) ), at = NULL, labels = paste( at ),
       type = "l", lwd = 2, lty = 1, col = par( "fg" ), log.ax = "y",
       las = 1, ann = FALSE, cex.ann = 0.8, xannx = 1/20,
       a.thin = seq( 1, length( age ), 2 ), ... )

Cplot( rates, age = as.numeric( rownames( rates ) ),
       per = as.numeric( colnames( rates ) ), grid = FALSE,
       c.grid = grid, ygrid = grid, col.grid = gray( 0.9 ),
       c.lim = NULL, ylim = range( rates[rates>0], na.rm=TRUE ),
       at = NULL, labels = paste( at ), c.lab = names( dimnames( rates ) )[2],
       ylab = deparse( substitute( rates ) ), type = "l", lwd = 2, lty = 1,
       col = par( "fg" ), log.ax = "y", las = 1, xannx = 1/20, ann = FALSE,
       cex.ann = 0.8, a.thin = seq( 1, length( age ), 2 ), ...  )

}
\arguments{
  \item{rates}{A two-dimensional table (or array) with rates to be plotted. It is
    assumed that the first dimension is age and the second is period.}
  \item{which}{A character vector with elements from
      \code{c("ap","ac","apc","pa","ca")}, indication which plots should
      be produced. One plot per element is produced. The first letter
      indicates the x-axis of the plot, the remaining which groups
      should be connected, i.e. \code{"pa"} will plot rates versus
      period and connect age-classes, and \code{"apc"} will plot rates
      versus age, and connect both periods and cohorts.}
  \item{age}{Numerical vector giving the means of the
    age-classes. Defaults to the rownames of \code{rates} as numeric.}
  \item{per}{Numerical vector giving the means of the periods. Defaults
    to the columnnames of \code{rates} as numeric.}
  \item{grid}{Logical indicating whether a background grid should be drawn.}
  \item{a.grid}{Logical indicating whether a background grid on the
      age-axis should be drawn. If numerical it indicates the
      age-coordinates of the grid.}
  \item{p.grid}{do. for the period.}
  \item{c.grid}{do. for the cohort.}
  \item{ygrid}{do. for the rate-dimension.}
  \item{col.grid}{The colour of the grid.}
  \item{a.lim}{Range for the age-axis.}
  \item{p.lim}{Range for the period-axis.}
  \item{c.lim}{Range for the cohort-axis.}
  \item{ylim}{Range for the y-axis (rates).}
  \item{at}{Position of labels on the y-axis (rates).}
  \item{labels}{Labels to put on the y-axis (rates).}
  \item{a.lab}{Text on the age-axis. Defaults to "Age".}
  \item{p.lab}{Text on the period-axis. Defaults to "Date of diagnosis".}
  \item{c.lab}{Text on the cohort-axis. Defaults to "Date of birth".}
  \item{ylab}{Text on the rate-axis. Defaults to the name of the rate-table.}
  \item{type}{How should the curves be plotted. Defaults to \code{"l"}.}
  \item{lwd}{Width of the lines. Defaults to 2.}
  \item{lty}{Which type of lines should be used. Defaults to 1, a solid line.}
  \item{log.ax}{Character with letters from \code{"apcyr"}, indicating
    which axes should be logarithmic. \code{"y"} and \code{"r"} both
    refer to the rate scale. Defaults to \code{"y"}.} 
  \item{las}{see \code{par}.}
  \item{ann}{Should the curves be annotated?}
  \item{a.ann}{Logical indicating whether age-curves should be annotated.}
  \item{p.ann}{do. for period-curves.}
  \item{c.ann}{do. for cohort-curves.}
  \item{xannx}{The fraction that the x-axis is expanded when curves are annotated.}
  \item{cex.ann}{Expansion factor for characters annotating curves.}
  \item{a.thin}{Vector of integers indicating which of the age-classes
    should be labelled.}
  \item{p.thin}{do. for the periods.}
  \item{c.thin}{do. for the cohorts.}
  \item{col}{Colours for the curves.}
  \item{a.col}{Colours for the age-curves.}
  \item{p.col}{do. for the period-curves.}
  \item{c.col}{do. for the cohort-curves.}
  \item{p.lines}{Should rates from the same period be connected?}
  \item{c.lines}{Should rates from the same cohort be connected?}
  \item{...}{Additional arguments pssed on to \code{matlines} when
    plotting the curves.}
}
\details{
Zero values of the rates are ignored. They are neiter in the plot nor in
the calculation of the axis ranges.
}  
\value{
  \code{NULL}. The function is used for its side-effect, the plot.
}
\author{
  Bendix Carstensen, Steno Diabetes Center,
  \url{http://bendixcarstensen.com}
}
\seealso{
  \code{\link{apc.frame}}
}
\examples{
data( blcaIT )
attach(blcaIT)

# Table of rates:
bl.rate <- tapply( D, list(age,period), sum ) /
           tapply( Y, list(age,period), sum )
bl.rate

# The four classical plots:
par( mfrow=c(2,2) )
rateplot( bl.rate*10^6 )

# The labels on the vertical axis could be nicer:
rateplot( bl.rate*10^6, at=10^(-1:3), labels=c(0.1,1,10,100,1000) ) 

# More bells an whistles
par( mfrow=c(1,3), mar=c(3,3,1,1), oma=c(0,3,0,0), mgp=c(3,1,0)/1.6 )
rateplot( bl.rate*10^6, ylab="", ann=TRUE, which=c("AC","PA","CA"),
                      at=10^(-1:3), labels=c(0.1,1,10,100,1000),
                      col=topo.colors(11), cex.ann=1.2 ) 
}
\keyword{hplot}
