
\name{sameLength}
\alias{sameLength}
\title{Make equal length of upstream and downstream regions in a profile list}
\description{
This function makes equal length of upstream and downstream regions in a profile list.
It is called by getProfiles if method ''middle'' is chosen and sameLength=T.
}
\usage{
sameLength(distribution, method=c("upstream", "downstream"))
}
\arguments{
\item{distribution}{a profile list}
\item{method}{character vector, indicating which regions of profiles should be made the same length. Only upstream and downstream are available.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
