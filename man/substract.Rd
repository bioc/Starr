
\name{substract}
\alias{substract}
\title{Substract mean or median from data}
\description{
This function substrcts the mean or median from each experiment. To use this normalization, the function normalize.Probes should be called with method="substract".
}
\usage{
substract(matrix, fun)
}
\arguments{
\item{matrix}{a matrix, containing the data}
\item{fun}{function, that should be applied to the experiments (mean or median)}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
