\name{lls}
\alias{lls}
\alias{clear}
\title{Functions to manage and explore the workspace
}
\description{These functions help you to find out what has gone wrong and to
start afresh if needed.
}
\usage{
lls(pos = 1, pat = "", all=FALSE, print=TRUE )
clear()
}
\arguments{
  \item{pos}{Numeric. What position in the search path do you want listed.}
  \item{pat}{Character. List only objects that have this string in their name.}
  \item{all}{Logical. Should invisible objects be printed too -
             see \code{\link{ls}} to which this argument is passed.}
  \item{print}{Logical. Should the result be printed?}
  }
\details{\code{lls} is designed to give a quick overview of the name, mode, class
  and dimension of the object in your workspace. They may not always be what you
  think they are.

  \code{clear} clears all your objects from workspace, and all attached objects
  too --- it only leaves the loaded packages in the search path; thus allowing a
  fresh start without closing and restarting R.
  }
\value{
\code{lls} returns a data frame with four character variables:
code{name},
code{mode},
code{class} and
code{size} and one row per object in the workspace (if \code{pos=1}).
\code{size} is either the length or the dimension of the object.
The data frame is by default printed with left-justified columns.
}
\author{\code{lls}: Unknown. Modified by Bendix Carstensen from a long
        forgotten snatch.
        
        \code{clear}: Michael Hills / David Clayton.}
\examples{
x <- 1:10
y <- rbinom(10, 1, 0.5)
m1 <- glm( y ~ x, family=binomial )
M <- matrix( 1:20, 4, 5 )
.M <- M
lls()
clear()
lls()
}
\keyword{attributes}