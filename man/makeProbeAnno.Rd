
\name{makeProbeAnno}
\alias{makeProbeAnno}
\title{ Creating a probeAnno object}
\description{
Creates a probeAnno object (package: Ringo) from a given Affymetrix bpmap file or a Nimblegen POS file. The posToProbeAnno function from the Ringo package is called to build the object.
}
\usage{
makeProbeAnno(posFile=NULL, bpmap=NULL, probeIDAsStrings=F)
}
\arguments{
\item{posFile}{path to the POS file}
\item{ bpmap}{ 	Either a list, created by the function readBpmap() from the affy package, or a path to the bpmap file. }
\item{probeIDAsStrings}{should the mapping of the probes to the rows in the assayData be integers or characters.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[Ringo:posToProbeAnnoEnvironment]{posToProbeAnno}},\code{\link[affxparser]{readBpmap}}}
\keyword{manip}
