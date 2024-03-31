\name{mh}
\alias{mh}
\alias{print.mh}
\title{
Mantel-Haenszel analyses of cohort and case-control studies
}
\description{
This function carries out Mantel-Haenszel comparisons in tabulated data 
derived from both cohort and case-control studies. 
}
\usage{
mh(cases, denom, compare=1, levels=c(1, 2), by=NULL,
     cohort=!is.integer(denom), confidence=0.9)
\method{print}{mh}(x, ...)
}
\arguments{
\item{cases}{
the table of case frequencies (a multiway array).
}
\item{denom}{
the denominator table. For cohort studies this should be a table of 
person-years observation, while for case-control studies it should be a 
table of control frequencies. 
}
\item{compare}{
the dimension of the table which defines the comparison groups (can be 
referred to either by number or by name). The default is the first 
dimension of the table.
}
\item{levels}{
a vector identifying (either by number or by name) the two groups to be 
compared. The default is the first two levels of the selected dimension.
}
\item{by}{
the dimensions not to be collapsed in the Mantel-Haenszel computations. 
Thus, this argument defines the structure of the resulting tables of
estimates and tests.  
}
\item{cohort}{
an indicator  whether the data derive from a cohort or a case-control 
study. If the denominator table is stored as an integer, a case-control
study is assumed.
}
\item{confidence}{
the approximate coverage probability for the confidence intervals to be 
computed.
}
\item{x}{a \code{mh} object
}
\item{...}{arguments passed on to \code{print}
}}
\value{
A list of class \code{mh} giving tables of rate (odds) ratio estimates,
their standard errors (on a log scale), lower and upper confidence
limits, chi-squared tests (1 degree of freedom) and the corresponding
p-values. The result list also includes numerator and denominator of the
Mantel-Haenszel estimates (q, r), and score test statistics and score
variance (u, v).
}
\section{Side Effects}{
None
}
\details{
Multiway tables of data
are accepted and any two levels of any dimension can be chosen as defining
the comparison groups. The rate (odds) ratio estimates and the associated
significance tests may be collapsed over all the remaining dimensions of the 
table, or over selected dimensions only, so that tables of estimates and 
tests are computed.
}
\references{
Clayton, D. and Hills, M. : Statistical Models in Epidemiology, Oxford 
University Press (1993). 
}
\seealso{
\code{\link{Lexis}}
}
\examples{
# If d and y are 3-way tables of cases and person-years 
# observation formed by tabulation by two confounders 
# (named "C1" and "C2") an exposure of interest ("E"), 
# the following command will calculate an overall 
# Mantel-Haenszel comparison of the first two exposure 
# groups.
#
# Generate some bogus data
dnam <- list( E=c("low","medium","high"), C1=letters[1:2], C2=LETTERS[1:4] )
d <- array( sample( 2:80, 24),
            dimnames=dnam, dim=sapply( dnam, length ) )
y <- array( abs( rnorm( 24, 227, 50 ) ),
            dimnames=dnam, dim=sapply( dnam, length ) )
mh(d, y, compare="E")
#
# Or, if exposure levels named "low" and "high" are to be 
# compared and these are not the first two levels of E :
#
mh(d, y, compare="E", levels=c("low", "high"))
#
# If we wish to carry out an analysis which controls for C1, 
# but examines the results at each level of C2:
#
mh(d, y, compare="E", by="C2")
#
# It is also possible to look at rate ratios for every 
# combination of C1 and C2 :
#
mh(d, y, compare="E", by=c("C1", "C2"))
#
# If dimensions and levels of the table are unnamed, they must 
# be referred to by number.
#
}
\keyword{htest}
