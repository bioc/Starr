
\name{plotGCbias}
\alias{plotGCbias}
\title{ Visualize GC-Bias of Hybridization }
\description{
Generates a plot showing the GC-bias of the hybridization.
}
\usage{
plotGCbias(intensity, sequence, main="")
}
\arguments{
\item{ intensity }{ a vector of type numeric, containing the measured intensities. }
\item{ sequence }{ a vector of type character, containing the sequences.}
\item{ main }{ head of the plot }
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[graphics]{boxplot}}}
\examples{
##
sequence <- unlist(lapply(1:50000, function(x) {paste(sample(c("A","T","C","G"),prob=c(0.3,0.3,0.2,0.2),25,replace=TRUE), collapse="")}))
values <- runif(50000,min=-2,max=2)
plotGCbias(values, sequence)
}
\keyword{hplot}