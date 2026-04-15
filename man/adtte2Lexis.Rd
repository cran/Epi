\name{adtte2Lexis {Epi}}
\alias{adtte2Lexis}

\title{
Create a Lexis Object Given an ADAM ADTTE Dataset.
}

\description{
Takes (a subset of) an ADTTE ADaM-format dataset along with choices
about recurrent events and outputs a corresponding
Lexis-object. Time since first entry to each state is
included as a time scale in the Lexis-object.
}

\usage{
adtte2Lexis(dat,
            start.state,
            Tr.ev,
            Term.st)
}

\arguments{
\item{dat}{ A data.frame with data on a format complying with CDISC
  ADaM BDS standards, see
  \url{https://www.cdisc.org/standards/foundational/adam/adam-basic-data-structure-bds-time-event-tte-analyses-v1-0}.
  See details for necessary columns, and Examples for an example.}
\item{start.state}{A character string, used as the name of the common starting
  state of all individuals.}
\item{Tr.ev}{A vector containing a subset of the values taken by the
  EVNTDSC variable in ADTTE indicating possible states.}
\item{Term.st}{A user-specified text string containing the name of the unique censoring state.}
}

\details{
The ADTTE data.frame, \code{dat}, should contain the columns
\code{USUBJID, STARTDT, ADT, AGE, EVNTDESC} and \code{CNSR}.  All other
columns, including baseline values and FL-variables, will be
discarded. If they are necessary in later analyses, they will need to be
added back to the resulting Lexis-object manually.
}

\value{
A Lexis-object with possible \code{lex.Cst} and \code{lex.Xst} states
specified in \code{Tr.ev} and \code{Tr.st}, and with time-scales for
calender period, \code{per}, age, \code{age}, time since entry,
\code{t.sE} and additional time-scales for time since first entry into
the states given in \code{Tr.ev}.  Baseline information is not retained.
}

\references{
\url{https://www.cdisc.org/standards/foundational/adam/adam-basic-data-structure-bds-time-event-tte-analyses-v1-0}
}

\author{
Lars Diaz, \email{larsjorgeful@gmail.com}
}

\examples{
## Creating columns in ADTTE dataset:
STUDYID <- c("ABC123", "ABC123", "ABC123", "ABC123",
             "ABC123", "ABC123", "ABC123", "ABC123", "ABC123", "ABC123",
             "ABC123") 

USUBJID <- c("SUBJ001", "SUBJ001", "SUBJ001", "SUBJ002",
             "SUBJ002", "SUBJ002", "SUBJ002", "SUBJ002", "SUBJ002",
             "SUBJ002", "SUBJ002") 

AGE <- c(65L, 65L, 65L, 67L, 67L, 67L, 67L, 67L, 67L, 67L,
         67L) 

SEX <- c("M", "M", "M", "F", "F", "F", "F", "F", "F", "F", "F") 

ARM <- c("Treatment A", "Treatment A", "Treatment A",
         "Treatment A", "Treatment A", "Treatment A", "Treatment A",
         "Treatment A", "Treatment A", "Treatment A", "Treatment A") 

PARAMCD <- c("HOSP", "RELAPSE", "RELAPSE", "HOSP", "HOSP",
             "HOSP", "HOSP", "RELAPSE", "RELAPSE", "RELAPSE", "RELAPSE") 

PARAM <- c("Time to Hospitl (days)",
           "Time to Relapse (days)",
           "Time to Relapse (days)",
           "Time to Hospitl (days)",
           "Time to Hospitl (days)",
           "Time to Hospitl (days)",
           "Time to Hospitl (days)",
           "Time to Relapse (days)",
           "Time to Relapse (days)",
           "Time to Relapse (days)",
           "Time to Relapse (days)") 

EVNTSEQ <- c(1, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4) 

STARTDT <- c("2023-01-01", "2023-01-01", "2023-01-31",
             "2023-01-03", "2023-02-13", "2023-05-10", "2023-06-29",
             "2023-01-03", "2023-03-18", "2023-04-26", "2023-06-18") 

ADT <- c("2023-07-20", "2023-01-31", "2023-07-20",
         "2023-02-13", "2023-05-10", "2023-06-29", "2023-07-22",
         "2023-03-18", "2023-04-26", "2023-06-18", "2023-07-22") 

ADTF <- c(NA_character_, NA_character_, NA_character_,
          NA_character_, NA_character_, NA_character_, NA_character_,
          NA_character_, NA_character_, NA_character_, NA_character_) 

AVAL <- c(200, 30, 170, 41, 86, 50, 23, 74, 39, 53, 34) 

CNSR <- c(1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1) 

EVNTDESC <- c("Censored", "Relapse", "Censored",
              "Hospitl", "Hospitl", "Hospitl",
              "Censored", "Relapse", "Relapse", "Relapse", "Censored") 

ANL01FL <- c("Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y", "Y",
             "Y", "Y") 

##collecting into example adtte dataset:
adtte <- data.frame(STUDYID, USUBJID, AGE, SEX, ARM, PARAMCD, PARAM, EVNTSEQ,
                    STARTDT, ADT, ADTF, AVAL, CNSR, EVNTDESC, ANL01FL)
head(adtte)

## transform to Lexis format
adtte2Lexis(adtte, 
            start.state = "Rand", 
            Tr.ev = c("Relapse", "Hospitl"),
            Term.st = "Censored") -> Lx
Lx
summary(Lx, timeScales = TRUE)
}
\keyword{manip}
