\name{diet}
\alias{diet}
\non_function{}
\title{Diet and heart data}
\description{
The \code{diet} data frame has 337 rows and 14 columns.
The data concern a subsample of subjects drawn from larger cohort
studies of the incidence of coronary heart disease (CHD). These subjects
had all completed a 7-day weighed dietary survey while taking part in
validation studies of dietary questionnaire methods. Upon the closure of
the MRC Social Medicine Unit, from where these studies were directed, it
was found that 46 CHD events had occurred in this group, thus allowing a
serendipitous study of the relationship between diet and the incidence
of CHD. 
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{id}{
subject identifier, a numeric vector
}
\item{doe}{
date of entry into follow-up study, a date
}
\item{dox}{
date of exit from the follow-up study, a date
}
\item{dob}{
date of birth, a date
}
\item{fail}{
status on exit, a numeric vector (codes 1, 3, 11, and 13 represent CHD events)
}
\item{job}{
occupation, a factor with levels
\code{Driver} 
\code{Conductor} 
\code{Bank worker} 
}
\item{month}{
month of dietary survey, a numeric vector
}
\item{energy}{
total energy intake (KCal per day/100), a numeric vector
}
\item{height}{
(cm), a numeric vector
}
\item{weight}{
(kg), a numeric vector
}
\item{fat}{
  fat intake (g/day), a numeric vector
}
\item{fibre}{
dietary fibre intake (g/day), a numeric vector
}
\item{energy.grp}{
high daily energy intake, a factor with levels
\code{<=2750 KCal}
\code{>2750 KCal}
}
\item{chd}{
CHD event, a numeric vector (1=CHD event, 0=no event)
}
}
\source{
The data are described and used extensively by Clayton and Hills,
Statistical Models in Epidemiology, Oxford University Press,
Oxford:1993. They were rescued from destruction by David Clayton and
reentered from paper printouts.
}
\examples{
data(diet)
}
\keyword{datasets}
