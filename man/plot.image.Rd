
\name{plotImage}
\alias{plotImage}
\title{ Reconstruct the array image}
\description{
Function to visualize spatial distribution of raw intensities on Affymetrix Oligoarrays.
}
\usage{
plotImage(cel)
}
\arguments{
\item{ cel }{ 	a character, specifying the path to the CEL file}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[affxparser]{readCel}},\code{\link[lattice]{levelplot}}}
\examples{
# dataPath <- system.file("extdata", package="Starr")
# plotImage(file.path(dataPath,"Rpb3_IP_chr1.cel"))
}
\keyword{hplot}