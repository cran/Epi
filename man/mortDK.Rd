\name{mortDK}
\alias{mortDK}
\docType{data}
\title{Population mortality rates for Denmark in 1-year age-classes.}
\description{
The \code{mortDK} data frame has 1820 rows and 21 columns.
}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{age}: \tab Age class, 0--89, 90:90+. \cr
\code{per}: \tab Calendar period, 38: 1938--42, 43: 1943--47, ..., 88:1988-92. \cr
\code{sex}: \tab Sex, 1: male, 2: female. \cr
\code{risk}: \tab Number of person-years in the Danish population. \cr
\code{dt}: \tab Number of deaths. \cr
\code{rt}: \tab Overall mortality rate in cases per 1000 person-years, i.e.
\code{rt=1000*dt/risk} \cr
           \tab Cause-specific mortality rates in cases per 1000 person-years: \cr
\code{r1}: \tab Infections \cr
\code{r2}: \tab Cancer. \cr
\code{r3}: \tab Tumors, benign, unspecific nature. \cr
\code{r4}: \tab Endocrine, metabolic. \cr
\code{r5}: \tab Blood. \cr
\code{r6}: \tab Nervous system, psychiatric. \cr
\code{r7}: \tab Cerebrovascular. \cr
\code{r8}: \tab Cardiac. \cr
\code{r9}: \tab Respiratory diseases, excl. cancer. \cr
\code{r10}: \tab Liver, excl. cancer. \cr
\code{r11}: \tab Digestive, other. \cr
\code{r12}: \tab Genitourinary. \cr
\code{r13}: \tab Ill-defined symptoms. \cr
\code{r14}: \tab All other, natural. \cr
\code{r15}: \tab Violent. \cr
}
}
\source{
Statistics Denmark, National board of health provided original data. 
Michael Andersson grouped the causes of death.
}
\examples{
data(mortDK)
}
\seealso{\code{\link{thoro}}, \code{\link{gmortDK}}}
\keyword{datasets}
