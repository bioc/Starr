

\name{getProfiles}
\alias{getProfiles}
\title{Get profiles of ChIP-signal over annotated features}
\description{
This function associates the measured ChIP signals to annotated features and stores the profile of each feature in a list.
Each profile is divided in three parts. The first entry is ''upstream'', which saves the signal upstream of start. Then
follows ''region'', which is from start to end and then ''downstream'', which stores the signals downstream of end.
}
\usage{
getProfiles(eSet, probeAnno, gffAnno, upstream, downstream, feature="ORF", borderNames, method, sameLength=T, fill=T, distance=8, spacing=4)
}
\arguments{
\item{eSet}{an ExpressionSet, containing on sample.}
\item{probeAnno}{a probeAnno object for the given ExpressionSet}
\item{gffAnno}{a data frame containing the annotation of the features of interest}
\item{upstream}{how many basepairs upstream of the feature start (feature start on the crick strand is end in gffAnno) should be taken.}
\item{downstream}{how many basepairs downtream of the feature start (feature end on the crick strand is start in gffAnno) should be taken.}
\item{feature}{name of the features (e.g. ORF, transcript, rRNA, ...)}
\item{borderNames}{names of the borders, flaking the feature (e.g. c("start", "stop"))}
\item{method}{Two methods are available. "middle", just takes the middle position of each probe and its corresponding value. This method should be used if the whole genome is tiled in an high resolution. "basewise" calculates for each base the mean of all probes overlapping with this position.}
\item{fill}{if "middle" is chosen the distance of the taken values equals the probe spacing on the chip. To avoid errors, because of regions lacking of probes, one can fill up these regions with NAs.}
\item{distance}{if method "middle" and fill==TRUE are chosen, distance is the max distance of no value occuring before filling in one NA.}
\item{spacing}{probe spacing on the chip. Only used for filling up with NAs in method "middle".}
\item{sameLength}{if method "middle" is chosen it can occur that the length of the upstream/downstream region vary a little. If sameLength==TRUE, upstream/downstream regions get all the same length.}
}
\value{
a list with the following entries
\item{ID}{the ID/name of the sample}
\item{upstream}{number of basepairs, taken upstream of the feature}
\item{downstream}{number of basepairs, taken upstream of the feature}
\item{method}{method used}
\item{borderNames}{names of the borders}
\item{feature}{feature type (e.g. "ORF")}
\item{profile}{a list which contains all profiles of the features in the gffAnno. Each entry consists of a list with the elements "upstream", "region", "downstream".}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link{fill}},\code{\link{fillNA}},\code{\link{mapFeatures}},\code{\link{getIntensities}},\code{\link{getFeature}},
	\code{\link{fill}},\code{\link[Starr]{getProfilesByBase}}}
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

# profile <- getProfiles(rpb3_rankpercentile_ratio, probeAnnoChr1, transcriptAnno, 500, 500, feature="transcript", borderNames=c("TSS", "TTS"), method="basewise", sameLength=T, fill=T, distance=8, spacing=4)

}
\keyword{ manip }
