\name{mod.Lexis}
\alias{glm.Lexis}
\alias{gam.Lexis}
\alias{coxph.Lexis}
\title{Fit intensity models to follow-up data in Lexis objects  
}
\description{
Modeling intensities based on Lexis objects, exploiting the structure of the
Lexis objects where the events and risk time have predefined
representations. This allows a simpler syntax than the
traditional explicit modeling using \code{\link{glm}}, \code{\link{gam}}
and \code{\link{coxph}}. But it is only a set of wrappers fro
\code{glm}, \code{gam} and \code{coxph}.
}
\usage{
  glm.Lexis( Lx,         # Lexis object	
           resp,         # 'to' states	
        formula,         # ~ model	
           xpos,         # 'from' states
           link = "log", # link function
          scale = 1,     # scaling of PY
        verbose = TRUE,  # report what is done?
             ... )
  gam.Lexis( Lx,         # Lexis object	
           resp,         # 'to' states	
        formula,         # ~ model	
           xpos,         # 'from' states
           link = "log", # link function
          scale = 1,     # scaling of PY
        verbose = TRUE,  # report what is done?
             ... )
coxph.Lexis( Lx, # Lexis object	
           resp, # 'to' states	
        formula, # timescale ~ model	
           xpos, # 'from' states
        verbose = TRUE,  # report what is done?
             ... )
}
\arguments{
  \item{Lx}{A \code{\link{Lexis}} object representing cohort follow-up.
}
  \item{resp}{Character vector of states \bold{to} which a transition is
    considered an event. May also be an integer vector in which case the
    reference will be to the position of levels of \code{lex.Xst}. 
}
  \item{formula}{Model formula describing the model for the
    intensity. For \code{glm} and \code{gam}, the formula should be
    one-sided; for \code{coxph} the formula should be two-sided and have
    the name of the time-scale used for baseline as the l.h.s.
}
\item{xpos}{Character vector of states \bold{from} which transitions
    are considered. May also be an integer vector in which case the
    reference will be to the position of levels of \code{lex.Cst}. If
    missing (that is not supplied), the entire \code{Lexis} object is
    used in the modeling.
}
\item{link}{Link function used, allowed values are \code{log} (the
    default), \code{identity} and \code{sqrt}, see the family
    \code{\link{poisreg}}. 
}
\item{scale}{Scalar. \code{lex.dur} is divided by this number before
  analysis, so that you can get resulting rates on a scale of your wish. 
}
  \item{verbose}{Print information on the states modeled?
}
  \item{\dots}{Arguments passed on to the methods.  }
}

\details{
  The \code{glm} and \code{gam} models are fitted using the family
  \code{\link{poisreg}} which is a bit faster than the traditional
  \code{poisson} family. The response variable for this family is a
  two-column vector of events and person-time respectively, so the
  predictions, for example using \code{\link{ci.pred}} does not require
  \code{lex.dur} (and would ignore this) as variable in the
  \code{newdata}. \code{ci.pred} returns the estimated rates in
  units of the \code{lex.dur} in the \code{Lexis} object, scaled by
  \code{scale}.

  Strictly speaking, it is a bit counter-intuitive to have the time-scale
  on the l.h.s. of the formula for the \code{coxph} since the time scale
  is also a predictor of the occurrence rate. On the other hand, calling
  \code{coxph} directly would also entail having the name of the time
  scale in the \code{Surv} object on the l.h.s. of the formula. So the
  inconsistency is merely carried over from \code{coxph}.

  The functions are slightly experimental so far. Argument names and
  ordering may change in the future. The \code{update} methods do not
  always work.
}

\value{\code{glm} returns a \code{\link{glm}} object, \code{gam} returns
  a \code{\link[mgcv]{gam}} object and \code{coxph} returns a
  \code{\link[survival]{coxph}} object. The returned objects all have an
  extra attribute, \code{Lexis}; a list with names \code{Exposure} and
  \code{Events}; character vectors of state names from which transitions
  are modeled and that are considered events, respectively. The
  \code{glm} object also has a \code{scale} element in the list, an scalar
  indicating the scaling of \code{lex.dur} before modeling. The
  \code{coxph} object also has a \code{Timescale} element in the list, a
  character indicating the underlying timescale variable.
}

\author{
Bendix Carstensen, \url{http://BendixCarstensen.com}.
}

\seealso{
\code{\link{Lexis}}, \code{\link{cutLexis}}
}

\examples{
library( Epi )
library( survival )
data( DMlate )

# Lexis object of total follow-up
mL <- Lexis( entry = list(age=dodm-dobth,per=dodm),
              exit = list(per=dox),
       exit.status = factor(!is.na(dodth),labels=c("Alive","Dead")),
              data = DMlate )

# Cut follow-up at start of insulin use
cL <- cutLexis( mL, cut = mL$doins,
              timescale = "per",
              new.state = "Ins",
       precursor.states = "Alive" )

# Split follow-up on age-axis
system.time( sL <- splitLexis( cL, breaks=0:25*4, time.scale="age") )
summary( sL )

# glm models for rates based on the time-split dataset by insulin and sex

# proportional hazards model with insulin as time-dependent variable
mt <- glm.Lexis( sL, resp="Dead",
                     ~ sex + lex.Cst + Ns(age,knots=c(15,3:8*10)) )

# prediction of mortality rates from "Alive" with and without PH assumption
nA <- data.frame( age=40:70, sex="M", lex.Cst="Alive" )
nI <- data.frame( age=40:70, sex="M", lex.Cst="Ins" )
matshade( nA$age, cbind( ci.pred(mt,nA),
                         ci.pred(mt,nI) )*1000, plot=TRUE,
          lwd=3, lty=1, log="y", col=c("black","blue","red"),
          xlab="Age", ylab="Mortality per 1000 PY" )
 
# gam models takes quite some time so we leave it out
\dontrun{
mt.gam <- gam.Lexis( sL, "Dead", ~ sex + lex.Cst + s(age), scale=1000 )
        }

# Fit a Cox model with age as baseline time scale and insulin as time-dependent
mt.cox <- coxph.Lexis( sL, "Dead", age ~ sex + lex.Cst )

# Pretty much the same results for regression paramters as the glm:
  ci.exp( mt    , subset="ex" )
# ci.exp( mt.gam, subset="ex" )
  ci.exp( mt.cox, subset="ex" )
}
\keyword{models}
