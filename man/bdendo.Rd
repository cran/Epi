\name{bdendo}
\alias{bdendo}
\non_function{}
\title{A case-control study of endometrial cancer}
\description{
The \code{bdendo} data frame has 315 rows and 13 columns.
These data concern a study in which each case of endometrial cancer was
matched with 4 controls. Matching was by date of birth (within one
year), marital status, and residence.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{set}{
Case-control set: a numeric vector
}
\item{d}{
Case or control: a numeric vector (1=case, 0=control)
}
\item{gall}{
Gall bladder disease: a factor with levels
\code{No} 
\code{Yes} 
}
\item{hyp}{
Hypertension: a factor with levels
\code{No} 
\code{Yes} 
}
\item{ob}{
Obesity: a factor with levels
\code{No} 
\code{Yes} 
}
\item{est}{
a factor with levels
\code{No} 
\code{Yes} 
}
\item{dur}{
Duration of conjugated oestrogen therapy: an ordered factor with levels
\code{0} < \code{1} < \code{2} < \code{3} < \code{4}
}
\item{non}{
Use of non oestrogen drugs: a factor with levels
\code{No} 
\code{Yes} 
}
\item{duration}{
Months of oestrogen therapy: a numeric vector
}
\item{age}{
a numeric vector
}
\item{cest}{
Conjugated oestrogen dose: an ordered factor with levels
\code{0} < \code{1} < \code{2} < \code{3}
}
\item{agegrp}{
a factor with levels
\code{55-59} 
\code{60-64} 
\code{65-69} 
\code{70-74} 
\code{75-79} 
\code{80-84} 
}
\item{age3}{
a factor with levels
\code{<64} 
\code{65-74} 
\code{75+} 
}
}
\source{
Breslow NE, and Day N, Statistical Methods in Cancer Research. Volume
I: The Analysis of Case-Control Studies. IARC Scientific
Publications, IARC:Lyon, 1980.  
}
\examples{
data(bdendo)
}
\keyword{datasets}
