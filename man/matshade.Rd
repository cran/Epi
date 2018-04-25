\name{matshade}
\alias{matshade}
\title{
  Plot confidence intervals as shaded areas
}
\description{
  Uses an x-vector and a matrix of 3*N columns with estimates and ci.s to
  produce the lines of estimates and shaded areas in transparent colours
  around the lines as confidence bands. 
}
\usage{
matshade(x, y, lty = 1, col = 1:5, alpha = 0.1, ...)
}
\arguments{
  \item{x}{Numerical vector. Unlike \code{matplot} this can only be a vector.
}
  \item{y}{A matrix with 3*N columns --- representing estimate and
    confidence bounds for N curves. Order of columns are assumed to be
    (est,lo,hi,est,lo,hi...) (or (est,hi,lo...))
}
  \item{lty}{Line types for the curves.
}
  \item{col}{Color of the estimated curves.
}
  \item{alpha}{Number in [0,1] indicating the transparency of the colors for
    the confidence intervals. Larger values makes the shades
    stronger. Can be a vector which then applies to the curves in turn.
}
  \item{\dots}{Arguments passed on to \code{\link{matlines}} for use
    when plotting the lines.
}
}
\details{All shaded areas are plotted first, the curves added
  afterwards. The function requires that a plotting region is active. 
}
\value{NULL. Used for its side effects.
}
\author{Bendix Carstensen, \url{BendixCarstensen.com}
}
\examples{
# Follow-up data of DM patients
data( DMlate )
mL <- Lexis( entry=list(age=dodm-dobth,per=dodm),
              exit=list(per=dox),
       exit.status=factor(!is.na(dodth),labels=c("Alive","Dead")),
              data=DMlate )
# Split follow-up and model by splines
sL <- splitLexis(mL,breaks=0:100)
mort <- glm( (lex.Xst=="Dead") ~ sex*Ns(age,knots=2:8*10),
             offset = log(lex.dur),
             family = poisson,
               data = sL )
# predict rates
ndM <- data.frame( age=10:90, sex="M", lex.dur=1000 )
ndF <- data.frame( age=10:90, sex="F", lex.dur=1000 )
prM <- ci.pred( mort, ndM )
prF <- ci.pred( mort, ndF )
# predict rate-ratio
MFr <- ci.exp( mort, ctr.mat=list(ndM,ndF) )
# empty plot frame
matplot( ndM$age, prM, col="transparent", log="y",
         xlab="Age", ylab="Mortality per 1000 PY (and RR)")
# put in lines and shades
matshade( ndM$age, cbind( MFr, prF, prM ), col=c(1,2,4), lwd=3 )
abline( h=1 )
}
\keyword{color}
