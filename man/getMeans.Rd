

\name{getMeans}
\alias{getMeans}
\title{Get mean ChIP-signal over annotated features}
\description{
\code{getMeans} calculates the mean ChIP-signal over annotated features
}
\usage{
getMeans(eSet, probeAnno, geneAnno, regions)
}
\arguments{
\item{eSet}{an ExpressionSet}
\item{probeAnno}{a probeAnno object for the given ExpressionSet}
\item{geneAnno}{a data frame containing the annotation of the features of interest}
\item{ regions }{ a data frame, containing four columns. The pos columns is a character with values out of c("start", "region", "end"). upstream and downstream ar integers, indicating how many bases upstream and downstream from the specified position in the feature are included. level is an integer, that says at which level the rectangle in the lower device should be plotted. The numeration goes from the bottom to the ceiling. cor is the correlation of the region, which is plotted in the upper panel.   }
}
\value{
a list. Each entry contains the mean signals over the specified regions (in the regions data frame) of all features in geneAnno. 
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link{getProfiles}}}
\examples{
## 
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

# transcriptAnno <- read.gffAnno(file.path(dataPath, "transcriptAnno.gff"), feature="transcript")
# filtered_orfs <- filterGenes(transcriptAnno, distance_us = 0, distance_ds = 0, minLength = 1000)

# pos <- c("start", "start", "start", "region", "region","region","region", "stop","stop","stop")
# upstream <- c(500, 0, 250, 0, 0, 500, 500, 500, 0, 250)
# downstream <- c(0, 500, 250, 0, 500, 0, 500, 0, 500, 250)
# info <- data.frame(pos=pos, upstream=upstream, downstream=downstream, stringsAsFactors=FALSE)
# means_rpb3 <- getMeans(rpb3_rankpercentile_ratio, probeAnnoChr1, transcriptAnno[which(transcriptAnno$name %in% filtered_orfs),], info)
}
\keyword{manip}
