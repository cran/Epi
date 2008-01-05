\name{S.typh}
\alias{S.typh}
\docType{data}
\title{Salmonella Typhimurium outbreak 1996 in Denmark.}
\description{
Matched case-control study of food poisoning.
}
\format{
  A data frame with 136 observations on the following 15 variables:
  \tabular{rl}{
    \code{id}: \tab Person identification \cr
    \code{set}: \tab Matched set indicator \cr
    \code{case}: \tab Case-control status (1:case, 0:control \cr
    \code{age}: \tab Age of individual \cr
    \code{sex}: \tab Sex of individual (1:male, 2:female) \cr 
    \code{abroad}: \tab Within the last two weeks visited abroad (1:yes, 0:no) \cr
    \code{beef}: \tab Within the last two weeks eaten beef \cr
    \code{pork}: \tab Within the last two weeks eaten pork \cr
    \code{veal}: \tab Within the last two weeks eaten veal \cr
    \code{poultry}: \tab Within the last two weeks eaten poultry \cr
    \code{liverp}: \tab Within the last two weeks eaten liverpaste \cr
    \code{veg}: \tab Within the last two weeks eaten vegetables \cr
    \code{fruit}: \tab Within the last two weeks eaten fruit \cr
    \code{egg}: \tab Within the last two weeks eaten eggs \cr
    \code{plant7}: \tab Within the last two weeks eaten meat from plant no. 7 \cr
  }
  }
\details{
  In the fall of 1996 an unusually large number of Salmonella
  Typhimurium cases were recorded in Fyn county in Denmark. The Danish
  Zoonosis Centre set up a matched case-control study to find the
  sources. Cases and two age-, sex- and residency-matched controls were
  telephone interviewed about their food intake during the last two
  weeks.

  The participants were asked at which retailer(s) they had purchased
  meat. Retailers were independently of this linked to meat processing
  plants, and thus participants were linked to meat processing
  plants. This way persons could be linked to (amongst other) plant no 7.} 
\source{
  Tine Hald.
}
\references{
  Molbak K and Hald T: Salmonella Typhimurium outbreak in late summer
  1996. A Case-control study. (In Danish:
  Salmonella typhimurium udbrud paa Fyn sensommeren
  1996. En case-kontrol undersogelse.) Ugeskrift for Laeger.,
  159(36):5372-7, 1997.
}
\examples{
data(S.typh)
}
\keyword{datasets}
