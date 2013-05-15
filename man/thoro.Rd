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
\describe{
\item{\code{id}}{Identification of person.}
\item{\code{sex}}{Sex, 1: male / 2: female.}
\item{\code{birthdat}}{Date of birth, \code{Date} variable.}
\item{\code{contrast}}{Group, 1: Thorotrast / 2: Control.}
\item{\code{injecdat}}{Date of contrast injection, \code{Date} variable.}
\item{\code{volume}}{Injected volume of Thorotrast in ml. Control
                  patients have a 0 in this variable.}
\item{\code{exitdat}}{Date of exit from the study, \code{Date} variable.}
\item{\code{exitstat}}{Status at exit, 1: dead / 2: alive,
                  censored at closing of study, 20 February 1992 / 3:
                  censored alive at some earlier date.}
\item{\code{cause}}{Cause of death. See causes in the helpfile for
                  \code{\link{gmortDK}}.}
\item{\code{liverdat}}{Date of liver cancer diagnosis, \code{Date} variable.}
\item{\code{liver}}{Indicator of liver cancer diagnosis. Not all livercancers
                  are histologically verified, hence
                  \code{liver >= hepcc + chola + hmang}}
\item{\code{hepcc}}{Hepatocellular carcinoma at \code{liverdat}.}
\item{\code{chola}}{Cholangiocellular carcinoma at \code{liverdat}.}
\item{\code{hmang}}{Haemangisarcoma carcinoma at \code{liverdat}.}
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
