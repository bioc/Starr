
\name{mapFeatures}
\alias{mapFeatures}
\title{Map middle positions of probes to annotated features}
\description{
This function maps the middle positions of probes to annotated features. It is called by the getProfiles method, if method
"middle" is chosen.
}
\usage{
mapFeatures(probeAnno, gffAnno, upstream, downstream, chr)
}
\arguments{
\item{probeAnno}{a probeAnno object for the given ExpressionSet}
\item{gffAnno}{a data frame containing the annotation of the features of interest}
\item{upstream}{how many basepairs upstream of the feature start (feature start on the crick strand is end in gffAnno) should be taken.}
\item{upstream}{how many basepairs downtream of the feature start (feature end on the crick strand is start in gffAnno) should be taken.}
\item{chr}{character vector of chromosomes, from which the features should be taken.}
}
\value{
a list with one entry for each feature, containing the mapping to probe indices in the assayData.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
