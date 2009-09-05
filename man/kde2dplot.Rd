
\name{kde2dplot}
\alias{kde2dplot}
\title{Compute density of a scatterplot}
\description{
A 2d density is computed by kde2D.
}
\usage{
kde2dplot(x,y,grid=50, ncol=30,nlevels=10,main="")
}
\arguments{
\item{x}{x coordinate of data }
\item{y}{y coordinate of data }
\item{grid}{Number of grid points in each direction}
\item{ncol}{number of colors}
\item{nlevels}{number of contour levels desired if levels is not supplied}
\item{main}{head of plot}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[MASS]{kde2d}}}
\examples{
##
points = 10^4
x <- rnorm(points/2)
x = c(x,x+2.5)
y <- x + rnorm(points,sd=0.8)
x = sign(x)*abs(x)^1.3
kde2dplot(x,y,main="Density Estimation")
}
\keyword{internal}