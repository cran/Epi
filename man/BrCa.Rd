\name{BrCa}
\alias{BrCa}
\docType{data}
\title{Clinical status,
relapse, metastasis and death in 2982 women with breast cancer.
}
\description{
This dataset is a transformation of the example dataset used by Crowther
and Lambert in their multistate paper.}
\usage{data(BrCa)}
\format{
  A data frame with 2982 observations on the following 17 variables:
  \describe{
    \item{\code{pid}}{Person-id; numeric}
    \item{\code{year}}{Calendar year of diagnosis}
    \item{\code{age}}{Age at diagnosis}
    \item{\code{meno}}{Menopausal status; a factor with levels \code{pre} \code{post}}
    \item{\code{size}}{Tumour size; a factor with levels \code{<=20 mm} \code{>20-50 mm} \code{>50 mm}}
    \item{\code{grade}}{Tumour grade; a factor with levels \code{2} \code{3}}
    \item{\code{nodes}}{Number of positive lymph nodes, a numeric vector}
    \item{\code{pr}}{Progesteron receptor level}
    \item{\code{pr.tr}}{Transformed progesteron level}
    \item{\code{er}}{Estrogen receptor level}
    \item{\code{hormon}}{Hormon therapy at diagnosis; a factor with levels \code{no} \code{yes}}
    \item{\code{chemo}}{Chemotherapy treatment; a factor with levels \code{no} \code{yes}}
    \item{\code{tor}}{Time of relapse, years since diagnosis}
    \item{\code{tom}}{Time of metastasis, years since diagnosis}
    \item{\code{tod}}{Time of death, years since diagnosis}
    \item{\code{tox}}{Time of exit from study, years since diagnosis}
    \item{\code{xst}}{Vital status at exit; a factor with levels \code{Alive} \code{Dead}}
  }
}
\details{
The dataset has been modified to contain the times (since diagnosis) of the events of
interest, to comply with the usual structure of data.
}
\source{
  The original data were extracted from:
  \url{http://fmwww.bc.edu/repec/bocode/m/multistate_example.dta}, this
  is modified representation of the same amount of information.
}
\references{
  The data were used as example in the paper by Crowther and
  Lambert: Parametric multistate survival models: Flexible
  modelling allowing transition-specific distributions with application
  to estimating clinically useful measures of effect differences; Stat
  Med 36 (29), pp 4719-4742, 2017. (No, it is not the paper, just the
  title.) 

  A parallel analysis using the \code{\link{Lexis}} machinery is available
  as: \url{http://bendixcarstensen.com/AdvCoh/papers/bcMS.pdf}
		}
\examples{
data(BrCa)
}
\keyword{datasets}
