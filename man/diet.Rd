\name{diet}
\alias{diet}
\docType{data}
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
\tabular{rl}{
\code{id}: \tab subject identifier, a numeric vector. \cr
\code{doe}: \tab date of entry into follow-up study, a
                 \code{\link{Date}} variable. \cr
\code{dox}: \tab date of exit from the follow-up study, a
                 \code{\link{Date}} variable. \cr
\code{dob}: \tab date of birth, a
                 \code{\link{Date}} variable. \cr
\code{y}: \tab number of years at risk, a numeric vector. \cr
\code{fail}: \tab status on exit, a numeric vector (codes 1, 3 and
13 represent CHD events) \cr
\code{job}: \tab occupation, a factor with levels
                 \code{Driver} 
                 \code{Conductor} 
                 \code{Bank worker} \cr
\code{month}: \tab month of dietary survey, a numeric vector \cr
\code{energy}: \tab total energy intake (kCal per day/100), a numeric
vector \cr
\code{height}: \tab (cm), a numeric vector \cr
\code{weight}: \tab (kg), a numeric vector \cr
\code{fat}: \tab   fat intake (10 g/day), a numeric vector \cr
\code{fibre}: \tab dietary fibre intake (10 g/day), a numeric vector \cr
\code{energy.grp}: \tab high daily energy intake, a factor with levels
\code{<=2750 KCal}
\code{>2750 KCal} \cr
\code{chd}: \tab CHD event, a numeric vector (1=CHD event, 0=no event) \cr
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
# Illustrate the follow-up in a Lexis diagram
Lexis.diagram( age=c(30,75), date=c(1965,1990),
               entry.date=cal.yr(doe), exit.date=cal.yr(dox), birth.date=cal.yr(dob), 
               fail=(fail>0), pch.fail=c(NA,16), col.fail=c(NA,"red"), cex.fail=1.0,
               data=diet )

}
\keyword{datasets}
