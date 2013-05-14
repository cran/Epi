\name{twoby2}
\alias{twoby2}
\title{Analysis of a two by two table}
\description{
  Computes the usual measures of association in a 2 by 2 table with
  confidence intervals. Also produces asymtotic and exact tests. Assumes
  that comparison of probability of the first column level between
  levels of the row variable is of interest. Output requires that the
  input matrix has meaningful row and column labels.
}
\usage{
twoby2(exposure, outcome,
       alpha = 0.05, print = TRUE, dec = 4,
       conf.level = 1-alpha, F.lim = 10000)
}
\arguments{
  \item{exposure}{If a table the analysis is based on the first two rows
    and first two columns of this. If a variable, this variable is
    tabulated against} 
  \item{outcome}{as the second variable}
  \item{alpha}{Significance level}
  \item{print}{Should the results be printed?}
  \item{dec}{Number of decimals in the printout.}
  \item{conf.level}{1-\code{alpha}}
  \item{F.lim}{If the table total exceeds \code{F.lim}, Fisher's exact
    test is not computed}
}
\value{A list with elements:
  \item{table}{The analysed 2 x 2 table augmented with probabilities and
    confidence intervals. The confidence intervals for the probabilities
    are computed using the normal approximation to the
    log-odds. Confidence intervals for the difference of proportions are
    computed using method 10 from Newcombe, Stat.Med. 1998, 17, pp.873
    ff.}
  \item{measures}{A table of Odds-ratios and relative risk with
    confidence intervals.}
  \item{p.value}{Exact p-value for the null hypothesis of OR=1}
}
\author{
  Mark Myatt. Modified by Bendix Carstensen.
}
\examples{
Treat <- sample(c("A","B"), 50, rep=TRUE )
Resp <- c("Yes","No")[1+rbinom(50,1,0.3+0.2*(Treat=="A"))]
twoby2( Treat, Resp )                 
twoby2( table( Treat, Resp )[,2:1] ) # Comparison the other way round
}
\keyword{univar}
\keyword{htest}
