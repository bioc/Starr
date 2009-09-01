
\name{getFeature}
\alias{getFeature}
\title{Getting profiles of one annotated features with a given mapping}
\description{
This function gets the probe intensities of a feature, that was mapped, using the mapFeatures function. 
It is called by the getIntensities function method, if method "middle" is chosen in getProfiles.
}
\usage{
getFeature(feature, k, mapping, expr_mat, gffAnno, reverse=T)
}
\arguments{
\item{feature}{The feature name}
\item{mapping}{output of the mapFeatures function}
\item{expr_mat}{the Expression matrix from the assayData}
\item{gffAnno}{a data frame, containing the gff annotation}
\item{k}{character or integer, indicating the chromosome}
\item{reverse}{if true, profiles on the crick strand are reversed, so they are in the right order}
}
\value{
a list, containing the profile for the feature.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
