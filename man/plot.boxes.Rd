
\name{plotBoxes}
\alias{plotBoxes}
\title{ boxplots of experiments }
\description{
Generates a boxplot of the of the given experiments.
}
\usage{
plotBoxes(eSet, col=NULL)
}
\arguments{
\item{ eSet }{Either an ExpressionSet or a matrix, containing the data.}
\item{col}{color, to fill the boxes}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[graphics]{boxplot}}}
\examples{
##
mat <- matrix(rnorm(1000000), ncol=2)
colnames(mat) <- c("Sample1", "Sample2")
mat[,1] <- mat[,1]-2
plotBoxes(mat)
}
\keyword{hplot}