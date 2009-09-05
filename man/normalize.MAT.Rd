

\name{normalize.MAT}
\alias{normalize.MAT}
\title{Normalize probes using MAT}
\description{
This function converts the ExpressionSet to an object of class tilingSet and calls the normalizeProbes function from the package rMAT. To use this normalization, the function normalize.Probes should be called with method="MAT".
}
\usage{
normalize.MAT(eSet, ...)
}
\arguments{
\item{eSet}{an ExpressionSet}
\item{...}{parameters, which are passed to the normalizeProbes function of the package rMAT}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[rMAT]{NormalizeProbes}}}
\keyword{internal}
