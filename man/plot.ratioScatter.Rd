
\name{plotRatioScatter}
\alias{plotRatioScatter}
\title{Plot ratios of all possible combinations of IP and CONTROL}
\description{
A matrix of pairwise scatterplots of the ratios is created. The lower panel shows the correlation of the data.
}
\usage{
plotRatioScatter(eSet, ip, control, density=F, sample=NULL, cluster=T)
}
\arguments{
\item{ eSet }{ 	an ExprssionSet or matrix, containing the data }
\item{ ip }{ an integer, or boolean vector, that indicates, which columns in the ExpressionSet are IP experiments  }
\item{ control }{ an integer, or boolean vector, that indicates, which columns in the ExpressionSet are CONTROL or REFERENCE experiments }
\item{ density }{ if TRUE, a density scatter plot is plotted. This plot shows the density of the data. }
\item{ sample }{ An integer, indicating the number of subsamples to take for the density scatterplot. This is only recommended if the data is very large, as the density computation takes some time. }#
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
z <- y*2
mat <- matrix(c(x,y,z), ncol=3)
colnames(mat) <- c("A", "B1", "B2")
plotRatioScatter(mat, c(TRUE, FALSE, FALSE), c(FALSE, TRUE, TRUE), density=TRUE)
}
\keyword{hplot}