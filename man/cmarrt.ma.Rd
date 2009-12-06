\name{cmarrt.ma}
\alias{cmarrt.ma}
\title{ Compute moving average statistics by incorporating the correlation structure}
\description{
 This function extends the moving average approach by incorporating the correlation structure. It also outputs the p-values of the standardized moving average statistics under the Gaussian approximation.
}
\usage{
cmarrt.ma(eSet, probeAnno, chr=NULL, M=NULL,frag.length,window.opt='fixed.probe')
}

\arguments{
  \item{eSet}{ExpressionSet containing the normalized ratio}
  \item{probeAnno}{probeAnno object with mapping}
  \item{chr}{which chromosome should be analysed? If chr==NULL, all chromosome in the probeAnno object are taken.}
  \item{M}{rough estimate of the percentage of bound probes. If unknown, leave it NULL.}
  \item{frag.length}{average fragment length from sonication.}
  \item{window.opt}{option for sliding window, either "fixed.probe" or "fixed.gen.dist". Default is 'fixed.probe'.}
}
\details{
Computation using \code{window.opt = "fixed.probe"} calculates the moving average statistics within a fixed number of probes and is more efficient. Use this option if the tiling array is regular with approximately constant resolution. \code{window.opt="fixed.gen.dist"} computes the moving average statistics over a fixed genomic distance. }
\value{
  \item{data.sort}{datafile sorted by genomic position.}
  \item{ma}{unstandardized moving average(MA) statistics.}
  \item{z.cmarrt}{standardized MA under correlation structure.}
  \item{z.indep}{standardized MA under independence (ignoring correlation structure).}
  \item{pv.cmarrt}{p-values of probes under correlation.}
  \item{pv.indep}{p-values of probes under independence (ignoring correlation structure).}
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }

\author{Pei Fen Kuan, Adam Hinz}
\note{The p-values are obtained under the Gaussian approximation. Therefore, it is important to check the normal quantile-quantile plot if the Gaussian approximation is valid. The function also outputs the computation under independence (ignoring the correlation structure) for comparisons.
}
\seealso{ \code{\link{plotcmarrt}},\code{\link{cmarrt.peak}} }
\examples{
# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))

# cels <- c(file.path(dataPath,"Rpb3_IP_chr1.cel"), file.path(dataPath,"wt_IP_chr1.cel"), 
# 	file.path(dataPath,"Rpb3_IP2_chr1.cel"))
# names <- c("rpb3_1", "wt_1","rpb3_2")
# type <- c("IP", "CONTROL", "IP")
# rpb3Chr1 <- readCelFile(bpmapChr1, cels, names, type, featureData=TRUE, log.it=TRUE)

# ips <- rpb3Chr1$type == "IP"
# controls <- rpb3Chr1$type == "CONTROL"

# rpb3_rankpercentile <- normalize.Probes(rpb3Chr1, method="rankpercentile")
# description <- c("Rpb3vsWT")
# rpb3_rankpercentile_ratio <- getRatio(rpb3_rankpercentile, ips, controls, description, fkt=median, featureData=FALSE)

# probeAnnoChr1 <- bpmapToProbeAnno(bpmapChr1)
# peaks <- cmarrt.ma(rpb3_rankpercentile_ratio, probeAnnoChr1, chr=NULL, M=NULL,250,window.opt='fixed.probe')
}
\keyword{ manip }
