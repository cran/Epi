\name{nickel}
\alias{nickel}
\docType{data}
\title{A Cohort of Nickel Smelters in South Wales}
\description{
The \code{nickel} data frame has 679 rows and 7 columns.
The data concern a cohort of nickel smelting workers in South Wales and
are taken from Breslow and Day, Volume 2. For comparison purposes,
England and Wales mortality rates (per 1,000,000 per annum)
from lung cancer (ICDs 162 and 163),
nasal cancer (ICD 160), and all causes, by age group and calendar period, are
supplied in the dataset \code{\link{ewrates}}.
}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{id}: \tab Subject identifier (numeric) \cr
\code{icd}: \tab ICD cause of death if dead, 0 otherwise (numeric) \cr
\code{exposure}: \tab Exposure index for workplace (numeric) \cr
\code{dob}: \tab Date of birth (numeric) \cr
\code{age1st}: \tab Age at first exposure (numeric) \cr
\code{agein}: \tab Age at start of follow-up (numeric) \cr
\code{ageout}: \tab Age at end of follow-up (numeric) \cr
}
}
\source{
Breslow NE, and Day N, Statistical Methods in Cancer Research. Volume
II: The Design and Analysis of Cohort Studies. IARC Scientific
Publications, IARC:Lyon, 1987.  
}
\examples{
data(nickel)
str(nickel)
}
\keyword{datasets}
