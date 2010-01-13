
\name{remap}
\alias{remap}
\title{Remap reporter sequences to the genome and create a new bpmap file}
\description{
This function remaps the reporter sequences on the chip on the genome and outputs a new bpmap annotation, containing only unique matches to the genome.
A remapping is recommended if the bpmap file was built on an outdated genome, or if sequences, that match the genome more than once should be excluded. 
}
\usage{
remap(bpmap=NULL, seqs=NULL, nseq=NULL, path="", complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE, return_bpmap=FALSE)
}
\arguments{
\item{ bpmap}{A list, created by the function readBpmap() from the affy package.}
\item{nseq}{Number of sequences, that are searched in one iteration.}
\item{seqs}{Sequences to search as a character vector}
\item{path}{path to genomic fasta files}
\item{complementary}{If TRUE, the sequences are searched in the complementary strand of the text}
\item{reverse}{If TRUE, the sequences are searched in the reverse strand of the text}
\item{reverse_complementary}{If TRUE, the sequences are searched in the reverse complementary strand of the text}
\item{return_bpmap}{If TRUE, the output is a list in bpmap format}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{
# dataPath <- system.file("extdata", package="Starr")

# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))
# newbpmap <- remap(bpmapChr1, nseq=5000000, path=dataPath, reverse_complementary=TRUE, return_bpmap=TRUE)
}
\keyword{manip}
