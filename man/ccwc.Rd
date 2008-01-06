\name{ccwc}
\alias{ccwc}
\title{Generate a nested case-control study}
\usage{
ccwc(entry=0, exit, fail, origin=0, controls=1, match=list(), include=list(), data=NULL, silent=F)
}
\arguments{
  \item{entry}{
    Time of entry to follow-up
  }
  \item{exit}{
    Time of exit from follow-up
  }
  \item{fail}{
    Status on exit (1=Fail, 0=Censored)
  }
  \item{origin}{
    Origin of analysis time scale
  }
  \item{controls}{
    The number of controls to be selected for each case
  }
  \item{match}{
    List of categorical variables on which to match cases and controls
  }
  \item{include}{
    List of other variables to be carried across into the case-control
    study
  }
  \item{data}{
    Data frame in which to look for input variables
  }
  \item{silent}{
    If  False, echos a . to the screen for each case-control set
    created; otherwise produces no output.
  }
}
\description{
Given the basic outcome variables for a cohort study: the time of entry 
to the cohort, the time of exit and the reason for exit ("failure" or
"censoring"), this function computes risk sets and generates a matched
case-control study in which each case is compared with a set of controls
randomly sampled from the appropriate risk set. Other variables may be
matched when selecting controls.
}
\value{
  The case-control study, as a  dataframe containing:
  \item{Set}{
    case-control set number
  }
  \item{Map}{
    row number of record in input dataframe
  }
  \item{Time}{
    failure time of the case in this set
  }
  \item{Fail}{
    failure status (1=case, 0=control)
  }
  These are followed by the matching variables, and finally by the
  variables in the \code{include} list
}
\references{
  Clayton and Hills, Statistical Models in Epidemiology, Oxford
  University Press, Oxford:1993.
  }
\author{
  David Clayton
  }
\seealso{
  \code{\link{Lexis}}
  }
\examples{
#
# For the diet and heart dataset, create a nested case-control study
# using the age scale and matching on job
#
data(diet)
dietcc <- ccwc(doe, dox, chd, origin=dob, controls=2, data=diet,
               include=energy, match=job)
}
\keyword{datagen}
