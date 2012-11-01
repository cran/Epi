\name{effx.match}
\alias{effx.match}
\title{Function to calculate effects for individually matched case-control studies}
\description{
The function calculates the effects of an exposure on a response,
possibly stratified by a stratifying variable, and/or controlled for one
or more confounding variables. 
}
\usage{
effx.match(response,
exposure,
match,
strata=NULL,
control=NULL,
base=1,
digits=3,
alpha=0.05,
data=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{response}{The \code{response} variable - must be numeric}
\item{exposure}{The \code{exposure} variable can be numeric or a factor}
\item{match}{The variable which identifies the matched sets}
\item{strata}{The \code{strata} stratifying variable - must be a factor}
\item{control}{ The \code{control} variable(s). These are passed as a
  list if there are more than one of them.} 
\item{base}{Baseline for the effects of a categorical exposure, default 1}
\item{digits}{Number of significant digits for the effects, default 3}
\item{alpha}{1 - confidence level}
\item{data}{\code{data} refers to the data used to evaluate the function}
}
\details{Effects are calculated odds ratios.
The function is a wrapper for clogit, from the survival package.
The k-1 effects for a categorical  exposure with k levels are relative 
to a baseline which, by default, is the first level. The effect of a metric (quantitative) 
exposure is calculated per unit of exposure.
The exposure variable can be numeric or a factor, but if it is an ordered factor the order will be ignored.
}
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
library(survival)
data(bdendo)

# d is the case-control variable, set is the matching variable.
# The variable est is a factor and refers to estrogen use (no,yes)
# The variable hyp is a factor with 2 levels and refers to hypertension (no, yes)
# effect of est on the odds of being a case
effx.match(d,exposure=est,match=set,data=bdendo)
# effect of est on the odds of being a case, stratified by hyp
effx.match(d,exposure=est,match=set,strata=hyp,data=bdendo)
# effect of est on the odds of being a case, controlled for hyp
effx.match(d,exposure=est,match=set,control=hyp,data=bdendo)
}
\keyword{ models }
\keyword{ regression }
