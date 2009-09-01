
\name{bpmapToProbeAnno}
\alias{bpmapToProbeAnno}
\title{ Creating a probeAnno object}
\description{
\code{bpmapToProbeAnno} This function allows the user to create a probeAnno environment that holds the mapping between
probes on the array and their genomic match position(s). As input, the function takes an Affymetrix bpmap file.
}
\usage{
bpmapToProbeAnno(bpmap, verbose=T, uniqueSeq=T)
}
\arguments{
\item{ bpmap}{ 	Either a list, created by the function readBpmap() from the affy package. Or a path to the bpmap file. }
\item{verbose}{ should the progress be printed out?}
\item{uniqueSeq}{ If TRUE, probes sequences that occur more than once on the chip (and consequently match several positions on the genome) are set to 1 in the probeAnno object. Unique probes are set to 0. If false, all probes are set to 0. To identify all unique and multiple matching probes, a remapping of the probes to the genome is recommended.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{
##
# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))

# probeAnnoChr1 <- bpmapToProbeAnno(bpmapChr1)
}
\keyword{manip}
