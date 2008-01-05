\name{bdendo}
\alias{bdendo}
\docType{data}
\title{A case-control study of endometrial cancer}
\description{
The \code{bdendo} data frame has 315 rows and 13 columns.
These data concern a study in which each case of endometrial cancer was
matched with 4 controls. Matching was by date of birth (within one
year), marital status, and residence.
}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{set}: \tab Case-control set: a numeric vector  \cr
\code{d}: \tab Case or control: a numeric vector (1=case, 0=control) \cr
\code{gall}: \tab Gall bladder disease: a factor with levels
                  \code{No} 
                  \code{Yes}. \cr
\code{hyp}: \tab Hypertension: a factor with levels
                 \code{No} 
                 \code{Yes}. \cr
\code{ob}: \tab Obesity: a factor with levels
                \code{No} 
                \code{Yes}. \cr
\code{est}: \tab A factor with levels
                 \code{No} 
                 \code{Yes}. \cr
\code{dur}: \tab Duration of conjugated oestrogen therapy: an ordered factor with levels
                 \code{0} < \code{1} < \code{2} < \code{3} < \code{4}. \cr
\code{non}: \tab Use of non oestrogen drugs: a factor with levels
                 \code{No} 
                 \code{Yes}. \cr
\code{duration}: \tab Months of oestrogen therapy: a numeric vector.  \cr
\code{age}: \tab A numeric vector. \cr
\code{cest}: \tab Conjugated oestrogen dose: an ordered factor with levels
\code{0} < \code{1} < \code{2} < \code{3}. \cr
\code{agegrp}: \tab A factor with levels
                    \code{55-59} 
                    \code{60-64} 
                    \code{65-69} 
                    \code{70-74} 
                    \code{75-79} 
                    \code{80-84} \cr
\code{age3}: \tab a factor with levels
                  \code{<64} 
                  \code{65-74} 
                  \code{75+}  \cr
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
