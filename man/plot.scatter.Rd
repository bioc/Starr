
\name{plotScatter}
\alias{plotScatter}
\title{ High level scatterplot of experiments}
\description{
A matrix of pairwise scatterplots is created. The lower panle shows the correlation of the data.
}
\usage{
plotScatter(eSet, density=F, cluster=T, sample=NULL)
}
\arguments{
\item{ eSet }{ 	an ExprssionSet or matrix, containing the data }
\item{ density }{ if TRUE, a density scatter plot is plotted. This plot shows the density of the data. }
\item{ sample }{ An integer, indicating the number of subsamples to take for the density scatterplot. This is only recommended if the data is very large, as the density computation takes some time. }
\item{cluster}{if cluster=T, the experiments are clustered and similiar experiments are plotted together.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[graphics]{pairs}}, \code{\link{densityscatter}}}
\examples{
##
points <- 10^4
x <- rnorm(points/2)
x <- c(x,x+2.5)
x <- sign(x)*abs(x)^1.3
y <- x + rnorm(points,sd=0.8)
mat <- matrix(c(x,y), ncol=2)
colnames(mat) <- c("a", "b")
plotScatter(mat, density=TRUE)
}
\keyword{hplot}