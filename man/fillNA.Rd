
\name{fillNA}
\alias{fillNA}
\title{Fill large spaces in profiles with NA}
\description{
This function fills spaces in the profiles, that exceed a given distance.
It is called by getProfiles, if method ''middle'' is chosen.
}
\usage{
fillNA(distribution, mapping, upstream, downstream, gffAnno, distance, spacing)
}
\arguments{
\item{distribution}{a list containing the profiles for each feature}
\item{mapping}{a mapping of, which was constructed by the mapFeatures function}
\item{gffAnno}{a data frame containing the annotation of the features of interest}
\item{upstream}{how many basepairs upstream of the feature start (feature start on the crick strand is end in gffAnno) should be taken.}
\item{downstream}{how many basepairs downtream of the feature start (feature end on the crick strand is start in gffAnno) should be taken.}
\item{distance}{spaces, exceeding this distance are filled with NA}
\item{spacing}{probe spacing on the chip}
}
\value{
a list with one entry for each feature.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
