
\name{getProfilesByBase}
\alias{getProfilesByBase}
\title{Get profiles of ChIP-signal over annotated features (basewise)}
\description{
This function associates the measured ChIP signals to annotated features and stores the profile of each feature in a list.
It is called by getProfiles, if method ''basewise'' is chosen. User should call this function.
}
\usage{
getProfilesByBase(eSet, probeAnno, chr, gffAnno, upstream, downstream)
}
\arguments{
\item{eSet}{an ExpressionSet, containing on sample.}
\item{probeAnno}{a probeAnno object for the given ExpressionSet}
\item{gffAnno}{a data frame containing the annotation of the features of interest}
\item{upstream}{how many basepairs upstream of the feature start (feature start on the crick strand is end in gffAnno) should be taken.}
\item{upstream}{how many basepairs downtream of the feature start (feature end on the crick strand is start in gffAnno) should be taken.}
\item{chr}{character vector of chromosomes, from which the features should be taken.}
}
\value{
a list with one entry for each feature.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
