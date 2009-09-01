
\name{read.gffAnno}
\alias{read.gffAnno}
\title{Reading gff annotation}
\description{
This functions reads the annotation from a gff file.
}
\usage{
read.gffAnno(gffFile, feature=NULL)
}
\arguments{
\item{gffFile}{path to file}
\item{feature}{feature to select ("character"). If feature="gene", then only rows, representing this feature are read.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{
##
# dataPath <- system.file("extdata", package="Starr")
# transcriptAnno <- read.gffAnno(file.path(dataPath, "transcriptAnno.gff"), feature="transcript")
}
\keyword{IO}
