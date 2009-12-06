\name{declare.bound}
\alias{declare.bound}
\title{Declare bound probes for a given error rate control}
\description{
  An internal function to be used by  \code{\link{cmarrt.peak}}.
}
\usage{
declare.bound(alpha, method, pv.cmarrt, pv.indep)
}
\arguments{
  \item{alpha}{error rate control for declaring bound region.}
  \item{method}{correction method inherited from \code{\link{p.adjust}}.}
  \item{pv.cmarrt}{a numeric vector of p-values for MA statistics obtained under correlation structure via Gaussian approximation.}
  \item{pv.indep}{a numeric vector of p-values for MA statistics obtained under independence via Gaussian approximation.}
}

\value{
   \item{bdd.cmarrt}{vector of 0's(unbound) or 1's(bound) for each probe under correlation.} 
   \item{bdd.indep}{vector of 0's(unbound) or 1's(bound) for each probe under independence.}
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }
\author{Pei Fen Kuan, Adam Hinz}
\seealso{ \code{\link{cmarrt.peak}} }
\keyword{internal}
