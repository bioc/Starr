\name{plotcmarrt}
\alias{plotcmarrt}
\title{ Histogram of p-values and normal QQ plots for standardized MA statistics}
\description{
 Plot the histograms of p-values and normal QQ plots under correlation structure and independence. 
}
\usage{
plotcmarrt(cmarrt.ma, freq=FALSE)
}
\arguments{
  \item{cmarrt.ma}{ output object from \code{\link{cmarrt.ma}}.}
  \item{freq}{see ?hist}
}
\details{
 Diagnostic plots for comparing the distribution of standardized MA statistics under correlation and independence. 
}
\value{
  Histogram of p-values and normal QQ plots under correlation structure and independence.
}
\references{P.F. Kuan, H. Chun, S. Keles (2008). CMARRT: A tool for the analysiz of ChIP-chip data from tiling arrays by incorporating the correlation structure. \emph{Pacific Symposium of Biocomputing}\bold{13}:515-526. }
\author{Pei Fen Kuan, Adam Hinz}
\note{If the normal quantile-quantile plot deviates from the reference line for unbound probes, this indicates that Gaussian approximation is not suitable for analyzing this data.}

\seealso{ \code{\link{cmarrt.ma}},\code{\link{qqnorm}} }
\examples{
# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))

# cels <- c(file.path(dataPath,"Rpb3_IP_chr1.cel"), file.path(dataPath,"wt_IP_chr1.cel"), 
 	file.path(dataPath,"Rpb3_IP2_chr1.cel"))
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

# plotcmarrt(peaks)
}
\keyword{hplot}