\name{nickel}
\alias{nickel}
\alias{ew.all}
\alias{ew.lung}
\alias{ew.nasal}
\non_function{}
\title{A Cohort of Nickel Smelters in South Wales}
\description{
The \code{nickel} data frame has 679 rows and 7 columns.
The data concern a cohort of nickel smelting workers in South Wales and
are taken from Breslow and Day, Volume 2. For comparison purposes,
England and Wales mortality rates (per 1,000,000 per annum)
from lung cancer (ICDs 162 and 163),
nasal cancer (ICD 160), and all causes, by age group and calendar period, are
supplied as the matrices \code{ew.lung}, \code{ew.nasal}, and
\code{ew.all} respectively. 
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{id}{
Subject identifier (numeric)
}
\item{icd}{
ICD cause of death if dead, 0 otherwise (numeric)
}
\item{exposure}{
Exposure index for workplace (numeric)
}
\item{dob}{
Date of birth (numeric)
}
\item{age1st}{
Age at first exposure (numeric)
}
\item{agein}{
Age at start of follow-up (numeric)
}
\item{ageout}{
Age at end of follow-up (numeric)
}
}
\source{
Breslow NE, and Day N, Statistical Methods in Cancer Research. Volume
II: The Design and Analysis of Cohort Studies. IARC Scientific
Publications, IARC:Lyon, 1987.  
}
\examples{
data(nickel)
}
\keyword{datasets}
