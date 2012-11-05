\name{brv}
\alias{brv}
\docType{data}
\title{Bereavement in an elderly cohort}
\description{
The \code{brv} data frame has 399 rows and 11 columns.
The data concern the possible effect of marital bereavement on
subsequent mortality. They arose from a survey of the physical and
mental health of a cohort of 75-year-olds in one large general
practice. These data concern mortality up to 1 January, 1990 (although
further follow-up has now taken place).

Subjects included all lived with a living spouse when they entered the
study. There are three distinct groups of such subjects: (1) those in
which both members of the couple were over 75 and therefore included in
the cohort, (2) those whose spouse was below 75 (and was not, therefore,
part of the main cohort study), and (3) those living in larger
households (that is, not just with their spouse).
}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{id}: \tab subject identifier, a numeric vector \cr
\code{couple}: \tab couple identifier, a numeric vector \cr
\code{dob}: \tab date of birth, a date \cr
\code{doe}: \tab date of entry into follow-up study, a date \cr
\code{dox}: \tab date of exit from follow-up study, a date \cr
\code{dosp}: \tab date of death of spouse, a date (if the spouse was still alive
                  at the end of follow-up, this was coded to January 1, 2000) \cr 
\code{fail}: \tab status at end of follow-up,
                  a numeric vector (0=alive,1=dead) \cr
\code{group}: \tab see Description, a numeric vector \cr
\code{disab}: \tab disability score, a numeric vector \cr
\code{health}: \tab perceived health status score, a numeric vector \cr
\code{sex}: \tab a factor with levels
                 \code{Male} 
                 \code{Female}
}
}
\source{
Jagger C, and Sutton CJ, Death after Marital Bereavement. Statistics in
Medicine, 10:395-404, 1991. (Data supplied by Carol Jagger).
}
\examples{
data(brv)
}
\keyword{datasets}
