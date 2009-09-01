
\name{plotPosBias}
\alias{plotPosBias}
\title{ Bias of hybridzation, depending on base position in sequence. }
\description{
\code{plotPosBias} generates a plot showing the bias of hybridzation, depending on base position in sequence.
}
\usage{
plotPosBias(intensity, sequence, main="", ylim)
}
\arguments{
\item{ intensity }{ a vector of type numeric, containing the measured intensities }
\item{ sequence }{ a vector of type character, containing the sequneces}
\item{ main }{ head of the plot }
\item{ylim}{ ylim of plot}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{
##
sequence <- unlist(lapply(1:50000, function(x) {paste(sample(c("A","T","C","G"),prob=c(0.3,0.3,0.2,0.2),25,replace=TRUE), collapse="")}))
values <- runif(50000,min=-2,max=2)
plotPosBias(values, sequence)
}
\keyword{hplot}