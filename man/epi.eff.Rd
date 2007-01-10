\name{epi.eff}
\alias{epi.eff}
\title{Function to calculate effects}
\description{
The function calculates the effects of an exposure on a response,
possibly stratified by a stratifying variable, and/or controlled for one
of more confounding variables.
}
\usage{
epi.eff(response, type = "metric", exposure, strata = NULL, control = NULL, fup = NULL, data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{The \code{response} variable - must be numeric}
  \item{type}{The type of response\code{type} - must be one of "metric",
  "binary", "failure", or "count"}
  \item{exposure}{The \code{exposure} variable can be numeric or a factor}
  \item{strata}{The \code{strata} stratifying variable - must be a factor}
  \item{control}{ The \code{control} variable(s) - must take the form
    data.frame(variable names)}
  \item{fup}{The \code{fup} variable contains the follow-up time for a
    failure response}. This must be numeric.
  \item{data}{\code{data} refers to the data frame used}
}
\details{
epi.eff(), without arguments, produces a list of arguments
}
\value{
%  ~Describe the value returned
 % If it is a LIST, use
  \item{comp1 }{Effects of exposure}
  \item{comp2 }{Tests of significance}
% ...
}
%\references{ ~put references to the literature/web site here ~ }
\author{Michael Hills}
%\note{ ~~further notes~~ }

 %~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
data(births)
attach(births)
epi.eff(response=bweight,type="metric",exposure=hyp)
epi.eff(response=lowbw,type="binary",exposure=hyp)
sex <- factor(sex)
epi.eff(response=bweight,type="metric",exposure=hyp,strata=sex)
epi.eff(response=lowbw,type="binary",exposure=hyp,strata=sex)
}
\keyword{models}
\keyword{regression}
