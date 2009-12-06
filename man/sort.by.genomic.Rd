\name{sortbygenomic}
\alias{sortbygenomic}
\title{Pre-process the data by genomic location}
\description{
 An internal function to be used by  \code{\link{cmarrt.ma}}. The input data is sorted by genomic location.
}
\usage{
sortbygenomic(chr, start, stop, logR)
}
\arguments{
  \item{chr}{ chromosome ID of each probe.}
  \item{start}{start coordinate of each probe.}
  \item{stop}{stop coordinate of each probe.}
  \item{logR}{probe measurements, e.g., log ratio.}
}
\value{
  \item{sort.data}{sorted array with 4 columns(regID,chr,start,stop,logR)}
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }
\author{Pei Fen Kuan, Adam Hinz}
\seealso{ \code{\link{cmarrt.ma}} }
\keyword{internal}