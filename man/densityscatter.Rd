
\name{densityscatter}
\alias{densityscatter}
\title{Compute density of a scatterplot}
\description{
A 2d density is computed by kde2D.
}
\usage{
densityscatter(x,y,pch=19,cex=0.1,ncol=30,grid=100,palette="heat", add=F,...)
}
\arguments{
\item{x}{x coordinate of data }
\item{y}{y coordinate of data }
\item{pch}{type of point}
\item{cex}{A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default}
\item{grid}{Number of grid points in each direction}
\item{ncol}{number of colors}
\item{palette}{color palette to choose}
\item{add}{should data points be added to an exisiting plot?}
\item{...}{parameters passed to plot or points}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[Starr]{kde2plot}},\code{\link[graphics]{plot.default}},
  \code{\link[graphics]{points.default}}}
\examples{
##
points = 10^4
x <- rnorm(points/2)
x = c(x,x+2.5)
y <- x + rnorm(points,sd=0.8)
x = sign(x)*abs(x)^1.3
densityscatter(x,y)
}
\keyword{hplot}