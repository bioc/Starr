
\name{writeWIG}
\alias{writeWIG}
\title{write ChIP-chip data to a *.wig file}
\description{
This function writes the all columns of the assayData to a wiggle file.
}
\usage{
writeWIG(expressionSet, probeAnno, file, chr=NULL, probeLength=NULL) 
}
\arguments{
\item{expressionSet}{an ExpressionSet object}
\item{probeAnno}{a probeAnno object}
\item{file}{path to write to}
\item{chr}{subset of chromosomes in probeAnno. If specified, only the subset is written to the file.}
\item{probeLength}{length of the probes on the chip.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{IO}
