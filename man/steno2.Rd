\name{steno2}
\alias{steno2}
\alias{st2clin}
\alias{st2alb}
\docType{data}
\title{
Clinical trial: Steno2 baseline and follow-up. 
}
\description{
Steno-2 was a clinical trial conducted at Steno Diabetes Center
1993-2001. The intervention was intensified treatment versus
conventional treatment of diabetes patients with micro-albuminuria. The
datsets here concern the extended follow-up of the trial population till
2015. Three files are provided: \code{steno2} with one record per
person, \code{st2clin} with one record per clinical visit and
\code{st2alb} with one record per transition between states of
albuminuria.

These dataset are entirely simulated, but designed to give
approximately the same results as the original.
}
\usage{data("steno2")
       data("st2clin")
       data("st2alb")
}
\format{

  \code{steno2} is a data frame with 160 observations on the following 14
  variables: 
  \describe{
    \item{\code{id}}{person id, numeric}
    \item{\code{allo}}{Original trial allocation, a factor with levels
      \code{Int} \code{Conv}} 
    \item{\code{sex}}{Sex, a factor with levels \code{F} \code{M}}
    \item{\code{baseCVD}}{0/1 indicator of preexisting CVD at baseline}
    \item{\code{deathCVD}}{0/1 indicator whether cause of death was CVD}
    \item{\code{doBth}}{Date of birth, a Date}
    \item{\code{doDM}}{Date of diabetes diagnosis, a Date}
    \item{\code{doBase}}{Date of entry to study, a Date}
    \item{\code{doCVD1}}{Date of 1st CVD event, a Date}
    \item{\code{doCVD2}}{Date of 2nd CVD event, a Date}
    \item{\code{doCVD3}}{Date of 3rd CVD event, a Date}
    \item{\code{doESRD}}{Date of end stage renal disease, a Date}
    \item{\code{doEnd}}{Date of exit from follow-up, a Date}
    \item{\code{doDth}}{Date of death, a Date}
  }

  \code{st2clin} is data frame with 750 observations on clinical
  measurements at different clinical visits:
  \describe{
    \item{\code{id}}{person id, numeric}
    \item{\code{doV}}{Date of clinical visit, a Date}
    \item{\code{a1c}}{Glycosylated hemoglobin, mmol/mol}
    \item{\code{chol}}{Total cholesterol, mg/mol}
    \item{\code{crea}}{Creatinine, mg/mol}
  }
  
  \code{st2alb} is data frame with 307 observations of changes in
  complication (albuminuria) state
  \describe{
    \item{\code{id}}{person id, numeric}
    \item{\code{doTr}}{Date of transition, a Date}
    \item{\code{state}}{State of albuminuria, factor with levels
                        \code{Norm}, \code{Mic}, \code{Mac}. All persons
                        begin in the state \code{Mic}ro-albuminuria.} 
  }

}
\details{
  The data are not the original; all values of measurements and dates
  have been randomly perturbed, to prevent identifiability of
  individuals. Analysis of these data will give only (very)
  approximately the same results as in the published article, and only
  some of the aspects of data are included.
}
\references{
  P. Gaede, J. Oellgaard, B. Carstensen, P. Rossing, H. Lund-Andersen,
  H. H. Parving & O. Pedersen: Years of life gained by multifactorial
  intervention in patients with type 2 diabetes mellitus and
  microalbuminuria: 21 years follow-up on the Steno-2
  randomised trial. Diabetologia (2016), 59, pp 2298-2307
}
\examples{
data(steno2)
data(st2alb)
L2 <- Lexis( entry = list(per = doBase,
                          age = doBase - doBth),
              exit = list(per = doEnd),
       exit.status = factor(deathCVD + !is.na(doDth),
                            labels=c("Mic","D(oth)","D(CVD)")),
                id = id,
              data = cal.yr(steno2) )
summary(L2)
#
# Cut at intermediate transitions
cut2 <- data.frame(lex.id = st2alb$id,
                      cut = cal.yr(st2alb$do),
                new.state = st2alb$state)
L3 <- rcutLexis(L2, cut2)
summary(L3)
#
# no direct transitions Mic <-> Mac allowed, so put a cut in between:
dd <- subset(L3, (lex.Cst == "Mac" & lex.Xst =="Norm") |
                 (lex.Cst =="Norm" & lex.Xst == "Mac"))
# artificial visits to the middle state Mic: 
cut3 <- data.frame( lex.id = dd$lex.id,
                       cut = dd$per + dd$lex.dur/2,
                 new.state = "Mic")
L4 <- rcutLexis(L3, cut3)
summary(L4)
#
# Show all transitions
boxes(L4, boxpos = list(x = c(15,15,15,85,85),
                        y = c(50,15,85,25,75)),
          show.BE = TRUE, scale.R = 1000,
          cex=0.8, pos.arr=0.7, font=1, font.arr=1)
}
\keyword{datasets}
