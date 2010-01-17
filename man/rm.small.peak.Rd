\name{rm.small.peak}
\alias{rm.small.peak}
\title{ Remove bound regions which consist of too few probes }
\description{
 An internal function to be used by  \code{\link{cmarrt.peak}}.
}
\usage{
rm.small.peak(bdd.method, minrun, chr, start, stop, pv)
}
\arguments{
  \item{bdd.method}{vector of 0's(unbound) or 1's(bound) for each probe. }
  \item{minrun}{minimum number of probes to be called a bound region.}
  \item{chr}{ chromosome ID of each probe.}
  \item{start}{start coordinate of each probe.}
  \item{stop}{stop coordinate of each probe.}
  \item{pv}{p-values for each probe.}	
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{bdd.method}{vector of 0's(unbound) or 1's(bound) for each probe.} 
 \item{size:}{number of probes within each bound region.}
 \item{peak.start}{first probe location for each bound region.}
 \item{peak.end}{last probe location for each bound region.}
 \item{chr.start}{start chromosome for each bound region.}
 \item{chr.stop}{stop chromosome for each bound region.}
 \item{loc.start}{start coordinate for each bound region.}
 \item{loc.stop}{stop coordinate for each bound region.}
 \item{minPv}{minimum p-value for each bound region.}
 \item{avePv}{average p-value for each bound region.}
	
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }

\author{Pei Fen Kuan, Adam Hinz}
\seealso{ \code{\link{cmarrt.peak}} }
\keyword{internal}
