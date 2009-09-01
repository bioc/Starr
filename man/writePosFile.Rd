
\name{writePosFile}
\alias{writePosFile}
\title{ Creating a pos file}
\description{
Writes a Nimblegen pos file from a given Affymetrix bpmap file.
}
\usage{
writePosFile(bpmap, file)
}
\arguments{
\item{ bpmap}{ 	Either a list, created by the function readBpmap() from the affy package. Or a path to the bpmap file. }
\item{ file }{ 	a character, specifying the path to the file to be written }
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{IO}
