

\name{expressionByFeature}
\alias{expressionByFeature}
\title{ Getting expression value by feature from an ExpressionSet}
\description{
This function gets the expression of a specified feature (e.g. orf, gene) from an ExpressionSet.
}
\usage{
expressionByFeature(eSet, fkt, method="median")
}
\arguments{
\item{ eSet}{ An ExpressionSet, containing the normalized expression values}
\item{fkt}{ Function to convert the featureNames (e.g. affy IDs) of eSet to the required features (e.g. ORFs)}
\item{method}{ If one feature (e.g. ORF) has more than one feature (e.g. affy ID) on the chip, the mean/median over the intensities is taken}
}
\seealso{\code{\link[base:get]{mget}}}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{manip}
