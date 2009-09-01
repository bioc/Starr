
\name{plotMA}
\alias{plotMA}
\title{ M versus A plot }
\description{
A matrix of M vs. A plots of each pair (ip, control) is produced.
}
\usage{
plotMA(eSet, ip=NULL, control=NULL, col=NULL)
}
\arguments{
\item{ eSet }{ 	an ExprssionSet or matrix, containing the data }
\item{ ip }{ an integer, or boolean vector, that indicates, which columns in the ExpressionSet are IP experiments  }
\item{ control }{ an integer, or boolean vector, that indicates, which columns in the ExpressionSet are CONTROL or REFERENCE experiments }
\item{col}{color, to fill the boxes}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[affy]{ma.plot}}}
\examples{
##
mat <- matrix(rnorm(1000000), ncol=4)
colnames(mat) <- c("Sample1", "Sample2", "Sample3", "Sample4")
mat[,1] <- mat[,1]^2
plotMA(mat, c(TRUE, FALSE, TRUE, FALSE), c(FALSE, TRUE, FALSE, TRUE))
}
\keyword{hplot}