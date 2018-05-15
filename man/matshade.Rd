\name{matshade}
\alias{matshade}
\title{
  Plot confidence intervals as shaded areas around lines.
}
\description{
  Uses an x-vector and a matrix of 3*N columns with estimates and ci.s
  to produce the lines of estimates and confidence intervals as shaded
  areas in transparent colours around the lines of the estimates.
}
\usage{
matshade( x, y, lty = 1,
          col = 1:(ncol(y)/3), col.shade=col, alpha=0.15,
         plot = dev.cur()==1,
          ... )
}
\arguments{
  \item{x}{Numerical vector. Unlike \code{\link{matplot}} this can only be a vector.
}
  \item{y}{A matrix with 3*N columns --- representing estimates and
    confidence bounds for N curves. Order of columns are assumed to be
    (est,lo,hi,est,lo,hi...) (or (est,hi,lo...))
}
  \item{lty}{Line types for the curves.
}
  \item{col}{Color(s) of the estimated curves.
}
  \item{col.shade}{Color(s) of the shaded areas. These are the colors
    that are made transparent by the \code{alpha} factor. Defaults to
    the same colors as the lines.
}
  \item{alpha}{Number in [0,1] indicating the transparency of the colors for
    the confidence intervals. Larger values makes the shades
    darker. Can be a vector which then applies to the curves in turn.
}
  \item{plot}{Logical. Should a new plot frame be started? If no device
    is active, the default is to start one, and plot all \code{y}s versus
    x in transparent color. On the rare occasion a device is open, but
    no plot have been called you will get an error telling that plot.new
    has not been called yet, in which case you should explicitly set
    \code{plot} to \code{TRUE}. 
}  
  \item{\dots}{Arguments passed on to \code{\link{matplot}} (if
    \code{plot=TRUE}) and \code{\link{matlines}} for use when plotting
    the lines. Note that \code{lwd=0} will cause lines to be omitted and
    only the shades be plotted. 
}
}

\details{All shaded areas are plotted first, the curves added
  afterwards, so that lines are not 'overshadowed'.

  Also not that if you repeat the same command, you will get the curves
  and the shaded areas overplotted in the same frame, so the effect is to
  have the shades darker, because the transparent colors are plotted on
  top of those from the first command.
}
\value{NULL. Used for its side effects.
}
\author{Bendix Carstensen, \url{BendixCarstensen.com}
}
\examples{
# Follow-up data of Danish DM patients
data( DMlate )
mL <- Lexis( entry=list(age=dodm-dobth,per=dodm),
              exit=list(per=dox),
       exit.status=factor(!is.na(dodth),labels=c("Alive","Dead")),
              data=DMlate )
# Split follow-up and model by splines
sL <- splitLexis( mL, breaks=0:100, time.scale="age")
\dontrun{
# the same thing with popEpi
sL <- splitMulti( mL, age=0:100 )    
}
mort <- glm( (lex.Xst=="Dead") ~ sex*Ns(age,knots=c(15,3:8*10)),
             offset = log(lex.dur),
             family = poisson,
               data = sL )
\dontrun{
# The counterpart with gam
library( mgcv )
mort <- gam( (lex.Xst=="Dead") ~ s(age,by=sex) + sex,
             offset = log(lex.dur),
             family = poisson,
               data = sL )
       }
# predict rates (per 1000 PY) for men and women
ndM <- data.frame( age=10:90, sex="M", lex.dur=1 )
ndF <- data.frame( age=10:90, sex="F", lex.dur=1 )
# gam objects ignores the offset in prediction
prM <- ci.pred( mort, ndM )*1000
prF <- ci.pred( mort, ndF )*1000
# predict rate-ratio
MFr <- ci.exp( mort, ctr.mat=list(ndM,ndF) )
# plot lines with shaded confidence limits
matshade( ndM$age, cbind( MFr, prF, prM ), col=c(1,2,4), lwd=3,
          log="y", xlab="Age", ylab="Mortality per 1000 PY (and RR)" )
abline( h=1 )
}
\keyword{color}
