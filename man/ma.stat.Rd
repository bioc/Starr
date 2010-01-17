\name{ma.stat}
\alias{ma.stat}
\title{Compute moving average statistics and p-values}
\description{
  An internal function to be used by  \code{\link{cmarrt.ma}}.
}
\usage{
ma.stat(regID, chr, start, stop, logR, frag.length, center, sigma, rho, window.opt = "fixed.probe")
}
\arguments{
  \item{regID}{consecutive genomic region ID.}
  \item{chr}{ chromosome ID of each probe.}
  \item{start}{start coordinate of each probe.}
  \item{stop}{stop coordinate of each probe.}
  \item{logR}{probe measurements, e.g., log ratio.}
 \item{frag.length}{average fragment length from sonication.}
  \item{center}{estimated center of the unbound probes from \code{\link{backgd.sd}}.}
  \item{sigma}{estimated standard deviation of the unbound probes \code{\link{backgd.sd}}.}
  \item{rho}{autocorrelation functions for maximum window size/lag of 30 probes \code{\link{backgd.sd}}.}
  \item{window.opt}{option for sliding window, either "fixed.probe" or "fixed.gen.dist". Default is 'fixed.probe'.}
}
\value{
  \item{ma}{unstandardized moving average(MA) statistics.}
  \item{z.cmarrt}{standardized MA under correlation structure.}
  \item{z.indep}{standardized MA under independence (ignoring correlation structure).}
  \item{pv.cmarrt}{p-values of probes under correlation.}
  \item{pv.indep}{p-values of probes under independence (ignoring correlation structure).}
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }
\author{Pei Fen Kuan, Adam Hinz}
\seealso{ \code{\link{cmarrt.ma}} }
\keyword{internal}