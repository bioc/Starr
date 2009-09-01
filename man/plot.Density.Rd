
\name{plotDensity}
\alias{plotDensity}
\title{ density plots of experiments }
\description{
Generates a plot, showing the densities of the experiments.
}
\usage{
plotDensity(eSet, oneDevice=T, main="")
}
\arguments{
\item{ eSet }{ 	an ExprssionSet or a matrix, containing the data }
\item{ oneDevice }{ should all lines be plotted to one device?}
\item{ main }{ head of the plot }
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[graphics]{plot.default}}, \code{\link[stats]{density}}}
\examples{
##
mat <- matrix(rnorm(1000000), ncol=2)
colnames(mat) <- c("Sample1", "Sample2")
mat[,1] <- mat[,1]-2
plotDensity(mat)
}
\keyword{hplot}