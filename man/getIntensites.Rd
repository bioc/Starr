
\name{getIntensities}
\alias{getIntensities}
\title{Get profile of anntated features from a probe mapping}
\description{
This function gets intensities from an ExpressionSet with a given mapping from the mapFeatures function.
It is called by getProfiles, if method ''middle'' is chosen.
}
\usage{
getIntensities(eSet, chr, mapping, gffAnno, upstream, downstream)
}
\arguments{
\item{eSet}{an ExpressionSet, containing on sample.}
\item{gffAnno}{a data frame containing the annotation of the features of interest}
\item{upstream}{how many basepairs upstream of the feature start (feature start on the crick strand is end in gffAnno) should be taken.}
\item{upstream}{how many basepairs downtream of the feature start (feature end on the crick strand is start in gffAnno) should be taken.}
\item{chr}{character vector of chromosomes, from which the features should be taken.}
\item{mapping}{a mapping of, which was constructed by the mapFeatures function}
}
\value{
a list with one entry for each feature.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
