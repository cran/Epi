\name{effx}
\alias{effx}
\title{Function to calculate effects}
\description{
The function calculates the effects of an exposure on a response,
possibly stratified by a stratifying variable, and/or controlled for one
or more confounding variables. 
}
\usage{
effx( response, type = "metric",
	         fup = NULL,     
	    exposure,     
	      strata = NULL,  
	     control = NULL,
             weights = NULL,
	       alpha = 0.05,   
	        base = 1,	      
	      digits = 3,     
	        data = NULL )    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{response}{The \code{response} variable - must be numeric}
\item{type}{The type of response\code{type} - must be one of "metric",
      "binary", "failure", or "count"}
\item{fup}{The \code{fup} variable contains the follow-up time for a
      failure response}. This must be numeric.  
\item{exposure}{The \code{exposure} variable can be numeric or a factor}
\item{strata}{The \code{strata} stratifying variable - must be a factor}
\item{control}{ The \code{control} variable(s) - these are passed as a
  list if there are more than one.} 
\item{weights}{Frequency weights for binary response only}
\item{base}{Baseline for the effects of a categorical exposure, default 1}
\item{digits}{Number of significant digits for the effects, default 3}
\item{alpha}{1 - confidence level}
\item{data}{\code{data} refers to the data used to evaluate the function}
}
\details{The function is a wrapper for glm. Effects are calculated as
  differences in means for a metric response, odds ratios for a binary
  response, and rate ratios for a failure or count response.  

  The k-1 effects for a categorical exposure with k levels are relative
  to a baseline which, by default, is the first level. The effect of a
  metric (quantitative) exposure is calculated per unit of exposure. 

  The exposure variable can be numeric or a factor, but if it is an
  ordered factor the order will be ignored.}
\value{
%  ~Describe the value returned
 % If it is a LIST, use
  \item{comp1 }{Effects of exposure}
  \item{comp2 }{Tests of significance}
% ...
}
\references{ www.mhills.pwp.blueyonder.co.uk }
\author{Michael Hills}
%\note{ ~~further notes~~ }

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
library(Epi)
data(births)
births$hyp <- factor(births$hyp,labels=c("normal","hyper"))
births$sex <- factor(births$sex,labels=c("M","F"))

# bweight is the birth weight of the baby in gms, and is a metric
# response (the default) 

# effect of hypertension on birth weight
effx(bweight,exposure=hyp,data=births) 
# effect of hypertension on birth weight stratified by sex
effx(bweight,exposure=hyp,strata=sex,data=births) 
# effect of hypertension on birth weight controlled for sex
effx(bweight,exposure=hyp,control=sex,data=births) 
# effect of gestation time on birth weight
effx(bweight,exposure=gestwks,data=births) 
# effect of gestation time on birth weight stratified by sex
effx(bweight,exposure=gestwks,strata=sex,data=births) 
# effect of gestation time on birth weight controlled for sex
effx(bweight,exposure=gestwks,control=sex,data=births) 

# lowbw is a binary response coded 1 for low birth weight and 0 otherwise
# effect of hypertension on low birth weight
effx(lowbw,type="binary",exposure=hyp,data=births)
# etc.
}
\keyword{ models }
\keyword{ regression }
