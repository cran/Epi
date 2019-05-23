\name{apc.LCa}
\alias{apc.LCa}
\alias{show.apc.LCa}
\title{Fit Age-Period-Cohort models and Lee-Carter models with effects
  modeled by natural splines.
}
\description{
\code{apc.LCa} fits an Age-Period-Cohort model and sub-models (using
\code{\link{apc.fit}}) as well as Lee-Carter models (using
\code{\link{LCa.fit}}). \code{show.apc.LCa} plots the models in little
boxes with their residual deviance with arrows showing their
relationships.
}
\usage{
apc.LCa( data,
  keep.models = FALSE,
          ... )
show.apc.LCa( x,
       dev.scale = TRUE,
             top = "Ad", ... )
}
\arguments{
  \item{data}{A data frame that must have columns \code{A}, \code{P},
    \code{D} and \code{Y}, see e.g. \code{\link{apc.fit}}
  }
  \item{keep.models}{Logical. Should the \code{apc} object and the 5
    \code{LCa} objects be returned too?
  }
  \item{...}{Further parameters passed on to \code{\link{LCa.fit}} or
    \code{\link{boxes.matrix}}. 
  }
  \item{x}{The result from a call to \code{apc.LCa}.}
  \item{dev.scale}{Should the vertical position of the boxes with the
    models be scales relative to the deviance between the Age-drift
    model and the extended Lee-Carter model?}
  \item{top}{The model presented at the top of the plot of boxes
    (together with any other model with larger deviance) when
    vertical position is scaled by deviances. Only "Ad", "AP", "AC",
    "APa" or "ACa" will make sense.}
}
\details{The function \code{apc.LCa} fits all 9 models (well, 10) available as
  extension and sub-models of the APC-model and compares them by
  returning deviance and residual df.
}
\value{A 9 by 2 matrix classified by model and deviance/df; optionally
  (if \code{models=TRUE}) a list with the matrix as \code{dev}, \code{apc}, an
  \code{apc} object (from \code{\link{apc.fit}}), and \code{LCa}, a list
  with 5 \code{LCa} objects (from \code{\link{LCa.fit}}). 
}
\author{
Bendix Carstensen, \url{http://BendixCarstensen.com}
}
\seealso{
\code{ \link{apc.fit}}, \code{\link{LCa.fit} }
}
\examples{
library( Epi )
clear()
# Danish lung cancer incidence in 5x5x5 Lexis triangles
data( lungDK )
lc <- subset( lungDK, Ax>40 )[,c("Ax","Px","D","Y")]
names( lc )[1:2] <- c("A","P")
head( lc )

al <- apc.LCa( lc, npar=c(9,6,6,6,10), keep.models=TRUE, maxit=500, eps=10e-3 )
show.apc.LCa( al, dev=TRUE )

# Danish mortality data
\dontrun{
data( M.dk )
mdk <- subset( M.dk, sex==1 )[,c("A","P","D","Y")]
head( mdk )

al <- apc.LCa( mdk, npar=c(15,15,20,6,6), maxit=50, eps=10e-3,
               quiet=FALSE, VC=FALSE )
show.apc.LCa( al, dev=FALSE )
show.apc.LCa( al, dev=TRUE )
show.apc.LCa( al, top="AP" )

# Fit a reasonable model to Danish mortality data and plot results
mAPa <- LCa.fit( mdk, model="APa", npar=c(15,15,20,6,6), c.ref=1930,
                 a.ref=70, quiet=FALSE, maxit=250 )
par( mfrow=c(1,3) )
plot( mAPa ) }
}
\keyword{regression}
\keyword{models}
