
\name{filterGenes}
\alias{filterGenes}
\title{ Filter Features/Genes }
\description{
This fucntion filters genes and other annotated features with respect to length, ovelaps and distance to other features.
}
\usage{
filterGenes(gffAnno, distance_us=500, distance_ds=500, minLength=-Inf, maxLength=Inf)
}
\arguments{
\item{gffAnno}{a data frame containing the annotation}
\item{distance_us}{how many basepairs upstream to the feature should not overlap with other features.}
\item{distance_ds}{how many basepairs downstream to the feature should not overlap with other features.}
\item{minLength}{minimal length of the feature}
\item{maxLength}{maximal length of the feature}
}
\value{
	a character vector with the names of the features, that passed the filter.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{
##
# dataPath <- system.file("extdata", package="Starr")
# transcriptAnno <- read.gffAnno(file.path(dataPath, "transcriptAnno.gff"), feature="transcript")
# filtered_transcripts <- filterGenes(transcriptAnno, distance_us = 0, distance_ds = 0, minLength = 1000)
}
\keyword{manip}
