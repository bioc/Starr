
\name{plotProfiles}
\alias{plotProfiles}
\title{Plotting ChIP profiles of one or more clusters}
\description{
\code{plotProfiles} plots the ChIP profiles of one or more clusters. Additionally on can display the distribution of e.g. gene expression in the clusters.
}
\usage{
plotProfiles(profiles, mfcol=NULL, mfrow=NULL, ylab="intensity", xlab="position", histograms=NULL, cluster, profileplot=T, meanprofile=T, ...)
}
\arguments{
\item{profiles}{a list constructed by the function getProfiles().}
\item{mfcol}{see ?par}
\item{mfrow}{see ?par}
\item{ylab}{see ?par}
\item{xlab}{see ?par}
\item{histograms}{a list of named vectors. Density plots are created for every vector and cluster.}
\item{cluster}{A named integer vector, that maps the features to the cluster.}
\item{profileplot}{should a clusterplot be shown?}
\item{meanprofile}{should the mean profiles of each cluster be plotted??}
\item{...}{arguments, passed to plot.default}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[stats]{density}}, \code{\link{profileplot}}}
\examples{
## 
sampls = 100
probes = 63
clus = matrix(rnorm(probes*sampls,sd=1),ncol=probes)
clus= rbind( t(t(clus)+sin(1:probes/10))+1:nrow(clus)/sampls , t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/sampls )
clustering = kmeans(clus,3)$cluster
names(clustering) <- 1:length(clustering)


profiles <- apply(clus, 1, function(x) {list(upstream=x[1:20], region=x[21:43], downstream=x[44:63])})
names(profiles) <- 1:length(clustering)
profiles <- list(profile=profiles, upstream=20, downstream=20, borderNames=c("start", "stop"))

plotProfiles(profiles, cluster=clustering, ylim=c(-1,2.5), type="l", lwd=2)
}
\keyword{hplot}
