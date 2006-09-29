\name{thoro}
\alias{thoro}
\docType{data}
\title{Thorotrast Study}
\description{
The \code{thoro} data frame has 2470 rows and 14 columns.
Each row represents one patient that have had cerebral angiography (X-ray of
the brain) with an injected contrast medium, either Thorotrast or
another one (the controls).
}
\format{
This data frame contains the following columns:
\tabular{rl}{
\code{id}:       \tab Identification of person. \cr
\code{sex}:      \tab Sex, 1: male / 2: female. \cr
\code{birthdat}: \tab Date of birth, \code{Date} variable. \cr
\code{contrast}: \tab Group, 1: Thorotrast / 2: Control. \cr
\code{injecdat}: \tab Date of contrast injection, \code{Date} variable. \cr
\code{volume}:   \tab Injected volume of Thorotrast in ml. Control
                  patients have a 0 in this variable. \cr
\code{exitdat}:  \tab Date of exit from the study, \code{Date} variable. \cr
\code{exitstat}: \tab Status at exit, 1: dead / 2: alive,
                  censored at closing of study, 20 February 1992 / 3:
                  censored alive at some earlier date. \cr
\code{cause}:    \tab Cause of death. See causes in the helpfile for
                  \code{\link{gmortDK}} \cr
\code{liverdat}: \tab Date of liver cancer diagnosis, \code{Date} variable. \cr
\code{liver}:    \tab Indicator of liver cancer diagnosis. Not all livercancers
                  are histologically verified, hence
                  \code{liver >= hepcc + chola + hmang} \cr
\code{hepcc}:    \tab Hepatocellular carcinoma at \code{liverdat}. \cr
\code{chola}:    \tab Cholangiocellular carcinoma at \code{liverdat}. \cr
\code{hmang}:    \tab Haemangisarcoma carcinoma at \code{liverdat}. \cr
}
}
\source{
M Andersson, M Vyberg, J Visfeldt, B Carstensen & HH Storm:
Primary liver tumours among Danish patients exposed to Thorotrast.
Radiation Research, 137, pp. 262--273, 1994.

M Andersson, B Carstensen HH Storm: Mortality and cancer
incidence after cerebral angiography. Radiation Research, 142,
pp. 305--320, 1995.
}
\examples{
data(thoro)
str(thoro)
}
\seealso{\code{\link{mortDK}}, \code{\link{gmortDK}}}
\keyword{datasets}
