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
and \code{\link{coxph}}. Requires that \code{lex.Cst} and \code{lex.Xst}
are defined as factors.

But it is just a set of wrappers fro
\code{glm}, \code{gam} and \code{coxph}.
}
\usage{
  glm.Lexis( Lx,         # Lexis object	
        formula,         # ~ model	
           from = preceding(Lx,to), # 'from' states
             to = absorbing(Lx)   , # 'to' states
         paired = FALSE, # only the pairwise
           link = "log", # link function
          scale = 1,     # scaling of PY
        verbose = TRUE,  # report what is done?
            \dots )        # further arguments to glm
  gam.Lexis( Lx,         # Lexis object	
        formula,         # ~ model	
           from = preceding(Lx,to), # 'from' states
             to = absorbing(Lx)   , # 'to' states	
         paired = FALSE, # only the pairwise
           link = "log", # link function
          scale = 1,     # scaling of PY
        verbose = TRUE,  # report what is done?
            \dots )        # further arguments to gam
coxph.Lexis( Lx,         # Lexis object	
        formula,         # timescale ~ model	
           from = preceding(Lx,to), # 'from' states
             to = absorbing(Lx)   , # 'to' states	
         paired = FALSE, # only the pairwise
        verbose = TRUE,  # report what is done?
            \dots )        # further arguments to coxph
}
\arguments{
  \item{Lx}{A \code{\link{Lexis}} object representing cohort follow-up.
}
  \item{formula}{Model formula describing the model for the
    intensity(-ies). For \code{glm} and \code{gam}, the formula should be
    one-sided; for \code{coxph} the formula should be two-sided and have
    the name of the time-scale used for baseline hazard as the l.h.s.
}
  \item{from}{Character vector of states \bold{from} which transitions
    are considered. May also be an integer vector in which case the
    reference will be to the position of levels of
    \code{lex.Cst}. Defaults to the collection of transient states
    immediately preceding the absorbing states. 
}
  \item{to}{Character vector of states \bold{to} which a transition is
    considered an event. May also be an integer vector in which case the
    reference will be to the position of levels of \code{lex.Xst}.
    Defaults to the set of absorbing states.  
}
  \item{paired}{Logical. Should the states mentioned in \code{to},
    rep. \code{from} be taken as pairs, indicating the only transitions
    modeled. If \code{FALSE} all transitions from any of the states in
    \code{from} to any states in \code{to} are modeled.
}
  \item{link}{Character; name of the link function used, allowed values
    are \code{'log'} (the default), \code{'identity'} and \code{'sqrt'},
    see the family \code{\link{poisreg}}. 
}
  \item{scale}{Scalar. \code{lex.dur} is divided by this number before 
    analysis, so that you can get resulting rates on a scale of your wish. 
}
  \item{verbose}{Print information on the states modeled?
}
  \item{\dots}{Further arguments passed on to \code{glm}, \code{glm} or
    \code{coxph} 
}
}

\details{
  The \code{glm} and \code{gam} models are fitted using the family
  \code{\link{poisreg}} which is a bit faster than the traditional
  \code{poisson} family. The response variable for this family is a
  two-column vector of events and person-time respectively, so the
  predictions, for example using \code{\link{ci.pred}} does not require
  \code{lex.dur} (and would ignore this) as variable in the
  \code{newdata}. \code{ci.pred} will return the estimated rates in
  units of the \code{lex.dur} in the \code{Lexis} object, scaled by
  \code{scale}, which has a default value of 1.

  The default is to model all transitions into any absorbing state by
  the same model (how wise is that??). If only \code{from} is given,
  \code{to} is set to all states reachable from \code{from}, which may
  be a really goofy model and if so a warning is issued. If only
  \code{to} is given, \code{from} is set to the collection of states
  from which \code{to} can be reached directly --- see
  \code{\link{preceding}} and its cousins. This convention means that if
  you have a \code{\link{Lexis}} object representing a simple survival
  analysis, with states, say, "alive" and "dead", you can dispense with
  the \code{from} and \code{to} arguments.

  Occasionally you only want to model a subset of the possible
  transitions from states in \code{from} to states in \code{to}, in
  which case you specify \code{from} and \code{to} as character vectors
  of the same length and set \code{paired=TRUE}. Then only transitions
  \code{from[i]} to \code{to[i]}, \code{i}=1,2,... will be modeled.
  
  There is no working \code{update} functions for these objects (yet). 
  
  Strictly speaking, it is a bit counter-intuitive to have the time-scale
  on the l.h.s. of the formula for the \code{coxph} since the time scale
  is also a predictor of the occurrence rate. On the other hand, calling
  \code{coxph} directly would also entail having the name of the time
  scale in the \code{Surv} object on the l.h.s. of the formula. So the
  inconsistency is merely carried over from \code{coxph}.
}

\value{\code{glm.Lexis} returns a \code{\link{glm}} object, which is
  also of class \code{glm.lex},
  \code{gam.Lexis} returns a \code{\link[mgcv]{gam}} object, which is
  also of class \code{gam.lex}, and
  \code{coxph.Lexis} returns a \code{\link[survival]{coxph}} object, which is
  also of class \code{coxph.lex}. These extra class attributes are meant
  to facilitate the (still pending) implementation of an \code{update} function.

  The returned objects all have an extra attribute, \code{Lexis} which
  is a list with entries
  \code{data}, the name of the \code{Lexis} object modeled (note that it
  is \emph{not} the object, only the name of it, which may not be portable);
  \code{trans}, a character vector of transitions modeled;
  \code{formula}, the model formula; and
  \code{scale}, the scaling applied to \code{lex.dur} before modeling. 
  
  Only the \code{glm} and \code{gam} objects have the \code{scale} element
  in the list; a scalar indicating the scaling of \code{lex.dur} before
  modeling. Note that the formula component of the \code{Lexis}
  attribute of a \code{coxph} object is a
  two-sided formula with the baseline time scale as the l.h.s.
}

\author{
Bendix Carstensen, \url{http://bendixcarstensen.com}.
}

\seealso{
  \code{\link{Lexis}},
  \code{\link{cutLexis}},
  \code{\link{mcutLexis}},
  \code{\link{addCov.Lexis}},
  \code{\link{absorbing}},
  \code{\link{transient}}
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
# ( consider splitMulti from the popEpi package )
summary( sL )

# glm models for rates based on the time-split dataset by insulin and sex

# Proportional hazards model with insulin as time-dependent variable
# - uses the defaul of modeling all transitions from both transient
# states ("Alive" and "Ins") to the absorbing state ("Dead"). 
mt <- glm.Lexis( sL, ~ sex + lex.Cst + Ns(age,knots=c(15,3:8*10)) )

# prediction of mortality rates from "Alive" with and without PH assumption
nA <- data.frame( age=40:70, sex="M", lex.Cst="Alive" )
nI <- data.frame( age=40:70, sex="M", lex.Cst="Ins" )
matshade( nA$age, cbind( ci.pred(mt,nA),
                         ci.pred(mt,nI) )*1000, plot=TRUE,
          lwd=3, lty=1, log="y", col=c("black","blue","red"),
          xlab="Age", ylab="Mortality per 1000 PY" )
 
# gam models may take some time to run so we leave it out
\dontrun{
mt.gam <- gam.Lexis( sL, ~ sex + lex.Cst + s(age), to="Dead",
                     scale=1000 )
        }

# Fit a Cox model for mortality with age as baseline time scale and
# insulin (lex.Cst) as time-dependent covariate 
mt.cox <- coxph.Lexis( sL, age ~ sex + lex.Cst, c("Alive","Ins"), "Dead" )

# Pretty much the same results for regression paramters as the glm:
  ci.exp( mt    , subset="ex" )
# ci.exp( mt.gam, subset="ex" )
  ci.exp( mt.cox, subset="ex" )
}
\keyword{models}
