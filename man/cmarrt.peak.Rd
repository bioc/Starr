\name{cmarrt.peak}
\alias{cmarrt.peak}
\title{ Obtain bound regions for a given error rate control}
\description{
  Obtain bound regions under a given error rate control using correction method from \code{\link{p.adjust}}. 
}
\usage{
cmarrt.peak(cmarrt.ma, alpha, method, minrun, asCherList=FALSE)
}
\arguments{
  \item{cmarrt.ma}{output object from \code{\link{cmarrt.ma}}.}
  \item{alpha}{error rate control for declaring bound region.}
  \item{method}{correction method inherited from \code{\link{p.adjust}}.}
  \item{minrun}{minimum number of probes to be called a bound region.}
  \item{asCherList}{If TRUE, result is returned as class cherList. See Ringo, for further description.}
}
\details{
  The function returns two objects, \code{cmarrt.bound} and \code{indep.bound}. Each object is a list of bound regions which can be accessed by \code{$chr} (chromosome), \code{$peak.start} (start coordinate of each bound region), \code{$peak.stop} (stop coordinate of each bound region),
 \code{$n.probe} (number of probes within each bound region), \code{$min.pv} (minimum p-values of each bound region), \code{$ave.pv} (average p-values of each bound region).}


\value{
  \item{cmarrt.bound}{list of bound regions obtained under correlation structure.}
  \item{indep.bound}{list of bound regions obtained under independence (ignoring correlation).}
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }
\author{Pei Fen Kuan, Adam Hinz}
\note{The list of bound regions obtained under independence (ignoring the correlation structure) is for comparison. It is not recommended to use this list for downstream analysis.}
\seealso{ \code{\link{cmarrt.ma}},\code{\link{p.adjust}} }
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
# peaklist <- cmarrt.peak(peaks)
}
\keyword{ manip }
