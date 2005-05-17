\name{lep}
\alias{lep}
\non_function{}
\title{An unmatched case-control study of leprosy incidence}
\description{
The \code{lep} data frame has 1370 rows and 7 columns. This was an
unmatched case-control study in which incident cases of leprosy in a
region of N. Malawi were compared with population controls. 
}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{id}: \tab subject identifier: a numeric vector \cr
\code{d}: \tab case/control status: a numeric vector (1=case, 0=control) \cr

\code{age}: \tab a factor with levels
\code{5-9} 
\code{10-14} 
\code{15-19} 
\code{20-24} 
\code{25-29} 
\code{30-44} 
\code{45+}  \cr

\code{sex}: \tab a factor with levels
\code{male},  
\code{female}  \cr 

\code{bcg}: \tab presence of vaccine scar, a factor with levels
\code{no} 
\code{yes}  \cr

\code{school}: \tab schooling, a factor with levels
\code{none} 
\code{1-5yrs} 
\code{6-8yrs} 
\code{sec/tert} \cr

\code{house}: \tab housing, a factor with levels
\code{brick} 
\code{sunbrick} 
\code{wattle} 
\code{temp}  \cr
}
}
\source{
The study is described in more detail in Clayton and Hills, Statistical
Models in Epidemiology, Oxford University Press, Oxford:1993.
}
\examples{
data(lep)
}
\keyword{datasets}
