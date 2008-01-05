\name{stat.table}
\alias{stat.table}
\alias{print.stat.table}
\title{Tables of summary statistics}
\description{
  \code{stat.table} creates tabular summaries of the data, using a
  limited set of functions. A list of index variables is used
  to cross-classify summary statistics.
}
\usage{
stat.table(index, contents = count(), data, margins = FALSE)
\method{print}{stat.table}(x, width=7, digits,...)
}
\arguments{
  \item{index}{A factor, or list of factors, used for cross-classification.
    If the list is named, then the names will be used when printing the
    table. This feature can be used to give informative labels to the
    variables.}
  \item{contents}{A function call, or list of function calls. Only a
    limited set of functions may be called (See Details below).  If the
    list is named, then the names will be used when printing the table.}
  \item{data}{an optional data frame containing the variables to be
    tabulated. If this is omitted, the variables will be searched for in the
    calling environment.}
  \item{margins}{a logical scalar or vector indicating which  marginal
    tables are to be calculated.  If a vector, it should be the same
    length as the \code{index} argument: values corresponding to
    \code{TRUE} will be retained in marginal tables.}
  \item{x}{an object of class \code{stat.table}.}
  \item{width}{a scalar giving the minimum column width when printing.}
  \item{digits}{a scalar, or named vector, giving the number
    of digits to print after the decimal point. If a named vector is used,
    the names should correspond to one of the permitted functions (See
    Details below) and all results obtained with that function will be
    printed with the same precision.}
  \item{...}{further arguments passed to other print methods.}
}
\details{
  This function is similar to \code{tapply}, with some enhancements:
  multiple summaries of multiple variables may be mixed in the
  same table; marginal tables may be calculated; columns and rows may
  be given informative labels; pretty printing may be controlled by the
  associated print method.

  This function is not a replacement for \code{tapply} as it also has
  some limitations.  The only functions that may be used in the
  \code{contents} argument  are: \code{\link{count}},
  \code{\link{mean}}, \code{\link{weighted.mean}}, \code{\link{sum}},
  \code{\link{quantile}}, \code{\link{median}}, \code{\link{IQR}},
  \code{\link{max}}, \code{\link{min}}, \code{\link{ratio}}, and
  \code{\link{percent}}.
  
  The \code{count()} function, which is the default, simply creates a
  contingency table of counts.  The other functions are applied to
  each cell created by combinations of the \code{index} variables.
}
\value{
  An object of class \code{stat.table}, which is a multi-dimensional
  array. A print method is available to create formatted one-way and
  two-way tables.
}
\author{Martyn Plummer}
\note{
  The permitted functions in the contents list %are overloaded functions that
  are defined inside \code{stat.table}.  They have the same interface as
  the functions callable from the command line, except for two
  differences. If there is an argument \code{na.rm} then its default
  value is always \code{TRUE}. A second difference is that the
  \code{quantile} function can only produce a single quantile in each call.
}

\seealso{\code{\link{table}}, \code{\link{tapply}},
\code{\link{mean}}, \code{\link{weighted.mean}},
\code{\link{sum}}, \code{\link{quantile}},
\code{\link{median}}, \code{\link{IQR}},
\code{\link{max}}, \code{\link{min}}, \code{\link{ratio}},
\code{\link{percent}}, \code{\link{count}}}
\examples{
data(warpbreaks)
# A one-way table
stat.table(tension,list(count(),mean(breaks)),data=warpbreaks)
# The same table with informative labels
stat.table(index=list("Tension level"=tension),list(N=count(),
           "mean number of breaks"=mean(breaks)),data=warpbreaks)

# A two-way table
stat.table(index=list(tension,wool),mean(breaks),data=warpbreaks)  
# The same table with margins over tension, but not wool
stat.table(index=list(tension,wool),mean(breaks),data=warpbreaks,
           margins=c(TRUE, FALSE))

# A table of column percentages
stat.table(list(tension,wool), percent(tension), data=warpbreaks)
# Cell percentages, with margins
stat.table(list(tension,wool),percent(tension,wool), margin=TRUE,
           data=warpbreaks)

# A table with multiple statistics
# Note how each statistic has its own default precision
a <- stat.table(index=list(wool,tension),
                contents=list(count(),mean(breaks),percent (wool)),
                data=warpbreaks)
print(a)
# Print the percentages rounded to the nearest integer
print(a, digits=c(percent=0))

# An Epidemiological example with follow-up time
data(nickel)
str(nickel)

# Make a grouped version of the exposure variable
nickel$egr <- cut( nickel$exposure, breaks=c(0, 0.5, 5, 10, 100), right=FALSE )
stat.table( egr, list( count(), percent(egr), mean( age1st ) ), data=nickel )

# Split the follow-up time by current age
nickel.ex <-
W.Lexis( entry=agein, exit=ageout, fail=icd \%in\% c(162,163),
         origin=0, breaks=seq(0,100,20),
         include=list( id, exposure, egr, age1st, icd ), data=nickel )
str( nickel.ex )

# Table of rates
stat.table( Time, list( n=count(), N=count(id), D=sum(Fail),
                        "Rate/1000"=ratio(Fail,Exit-Entry,1000) ),
            margin=1, data=nickel.ex )
# Two-way table of rates and no. persons contributing
stat.table( list(age=Time, Exposure=egr),
            list( N=count(id), D=sum(Fail), Y=sum((Exit-Entry)/1000),
                  Rate=ratio(Fail,Exit-Entry,1000) ),
            margin=TRUE, data=nickel.ex )
}
\keyword{iteration}
\keyword{category}
