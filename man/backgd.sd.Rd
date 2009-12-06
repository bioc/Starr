\name{backgd.sd}
\alias{backgd.sd}
\title{Background parameters (internal function)}
\description{
  An internal function to be used by  \code{\link{cmarrt.ma}}.
}
\usage{
backgd.sd(regID, logR, M = NULL, chr, start, stop)
}
\arguments{
  \item{regID}{consecutive genomic region ID.}
  \item{logR}{probe measurements, e.g., log ratios, t-statistics.}
  \item{M}{rough estimate of the percentage of bound probes. If unknown, leave it NULL.}
\item{chr}{ chromosome ID of each probe.}
  \item{start}{start coordinate of each probe.}
  \item{stop}{stop coordinate of each probe.}	
}

\value{
  \item{center}{estimated center of the unbound probes.}
  \item{sigma}{estimated standard deviation of the unbound probes.}
  \item{rho}{autocorrelation functions for maximum window size/lag of 30 probes.}
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }
\author{Pei Fen Kuan, Adam Hinz}
\seealso{ \code{\link{cmarrt.ma}} }
\keyword{internal}
