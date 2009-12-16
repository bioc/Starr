
\name{remap}
\alias{remap}
\title{Remap reporter sequences to the genome and create a new bpmap file}
\description{
This function remaps the reporter sequences on the chip on the genome and outputs a new bpmap annotation, containing only unique matches to the genome.
A remapping is recommended if the bpmap file was built on an outdated genome, or if sequences, that match the genome more than once should be excluded. 
}
\usage{
remap(bpmap, nseq, bsg, chromosome_wise=FALSE, mySeqs=NULL, complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE)
}
\arguments{
\item{ bpmap}{A list, created by the function readBpmap() from the affy package.}
\item{nseq}{Number of sequences, that are searched in one iteration.}
\item{chromosome_wise}{If the genome very large (e.g. mouse or human), the sequences should be searched chromosome_wise and not at once in all chromosomes to safe memory.}
\item{bsg}{An object of class BSGenome, which contains the genomic sequences.}
\item{mySeqs}{Names of the genomic sequences in the BSGenome object, to analyse. If NULL, all sequences are taken.}
\item{complementary}{If TRUE, the sequences are searched in the complementary strand of the text}
\item{reverse}{If TRUE, the sequences are searched in the reverse strand of the text}
\item{reverse_complementary}{If TRUE, the sequences are searched in the reverse complementary strand of the text}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{

# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath,"Scerevisiae_tlg_chr1.bpmap"))

# newbpmap <- remap(bpmapChr1, 500000, Scerevisiae)

## Create a binary bpmap file using functions from the affy package ...
# writeTpmap("Scerevisiae_new.tpmap", newbpmap)
# tpmap2bpmap("Scerevisiae_new.tpmap", "Scerevisiae_new.bpmap")

## ... or create a probeAnno object
# pA <- bpmapToProbeAnno(newbpmap)
}
\keyword{manip}
