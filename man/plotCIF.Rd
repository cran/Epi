\name{plotCIF}
\alias{plotCIF}
\alias{stackedCIF}

\title{Plotting Aalen-Johansen curves for competing events
}

\description{Function \code{plotCIF} plots, for one or more groups, the
cumulative incidence curves for a selected event out of two or more
competing events. Function \code{stackedCIF} plots, for one group or
population, the cumulative incidence curves for two or more competing
events such that the cumulative incidences are stacked upon each
other. The CIFs are are estimated by the Aalen-Johansen method.
}
\usage{
## S3 method for class 'survfit'
 plotCIF( x, event = 1,
              xlab = "Time",
              ylab = "Cumulative incidence",
              ylim = c(0, 1),
               lty = 1,
               col = "black", ... )

## S3 method for class 'survfit'
stackedCIF( x, group = 1,
                 col = "black",
                fill = "white",
                ylim = c(0,1),
                xlab = "Time",
                ylab = "Cumulative incidence", ... )
}
\arguments{
  
\item{x}{An object of class \code{\link[survival]{survfit}}, the \code{type} of
  \code{event} in \code{Surv()} being "\code{mstate}"; the first level
  of the event factor represents censoring and the remaining ones the
  alternative competing events. 
}
\item{event}{Determines the event for which the cumulative incidence
  curve is plotted by \code{plotCIF}.
}
\item{group}{An integer showing the selected level of a possible
  grouping factor appearing in the model formula in \code{survfit} when
  plotting by \code{stackedCIF} 
}
\item{col}{A vector specifying the plotting color(s) of the curve(s) for
  the different groups in \code{\link{plotCIF}}-- default: all "black". 
}
\item{fill}{A vector indicating the colours to be used for shading the
  areas pertinent to the separate outcomes in \code{\link{stackedCIF}}
  - default: all \code{"white"}.
}
\item{xlab}{Label for the $x$-axis.
}
\item{ylab}{Label for the $y$-axis.
}
\item{ylim}{Limits of the $y$-axis.
}
\item{lty}{A vector specifying the line type(s) of the curve(s) for
    the different groups - default: all 1 (=solid). 
}
\item{\dots}{Further graphical parameters to be passed.
}
}
\details{
The order in which the curves with \code{\link{stackedCIF}} are piled
upon each other is the same as the ordering of the values or levels of
the competing events in the pertinent event variable. The ordering can
be changed by permuting the levels as desired using function
\code{Relevel}, after which \code{survfit} is called with the relevelled
\code{event} variable in \code{Surv()}
}
\value{No value is returned but a plot is produced as a side-effect.
}
\references{Putter, H., Fiocco, M., Geskus, R.B. (2007). 
Tutorial in biostatistics: competing risks and multi-state models. 
Statistics in Medicine, 26: 2389--2430.
}
\author{Esa Laara, \email{esa.laara@oulu.fi}
}
\note{
Aalen-Johansen curves for competing events in several groups can also
be plotted by function \code{\link[survival]{plot.survfit}} of the survival
library as well as by some functions in other packages covering analysis
of time-to-event data.}
\seealso{
  \code{\link[survival]{survfit}},
  \code{\link{plot}},
  \code{\link[survival]{plot.survfit}}. 
}
\examples{
library(survival)   #  requires version 2.39-4 or later
head(mgus1)
#  Aalen-Johansen estimates of CIF are plotted by sex for two 
#  competing events: (1) progression (pcm), and (2) death, in 
#  a cohort of patients with monoclonal gammopathy.

#  The data are actually covering transitions from pcm to death, too,
#  for those entering the state of pcm. Such patients have two rows
#  in the data frame, and in their 2nd row the 'start' time is 
#  the time to pcm (in days). 

#  In our analysis we shall only include those time intervals with value 0
#  for variable 'start'. Thus, the relevant follow-up time is represented 
#  by variable 'stop' (days). For convenience, days are converted to years.

fitCI <- survfit(Surv(stop/365.25, event, type="mstate") ~ sex,
              data= subset(mgus1, start==0) )
par(mfrow=c(1,2))
plotCIF(fitCI, event = 1, col = c("red", "blue"),
  main = "Progression", xlab="Time (years)" )
text( 38, 0.15, "Men", pos = 2)
text( 38, 0.4, "Women", pos = 2)
plotCIF(fitCI, event = 2, col = c("red", "blue"), 
  main = "Death", xlab="Time (years)" )
text( 38, 0.8, "Men", pos = 2)
text( 38, 0.5, "Women", pos = 2)

par(mfrow=c(1,2))
stackedCIF(fitCI, group = 1, fill = c("gray80", "gray90"),
  main = "Women", xlab="Time (years)" )	
text( 36, 0.15, "PCM", pos = 2)
text( 36, 0.6, "Death", pos = 2)
stackedCIF(fitCI, group = 2, fill = c("gray80", "gray90"), 
  main = "Men", xlab="Time (years)" )
text( 39, 0.10, "PCM", pos = 2)
text( 39, 0.6, "Death", pos = 2)	
}
