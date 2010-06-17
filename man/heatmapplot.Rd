
\name{heatmapplot}
\alias{heatmapplot}
\title{heatmapplot}
\description{
Heatmap representation of binding profiles
}
\usage{
heatmapplot(profiles, colpal=c("black","dark blue","dark green", "green","gold", "yellow"), abl=NULL, subset=NULL)
}
\arguments{
\item{profiles}{a list of profiles returned by getProfiles(). Features must have same lengths.}
\item{colpal}{color palette for intensity coding}
\item{abl}{positions of vertical lines that are added to the panel}
\item{subset}{subset of genes in the list that should be plotted}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{hplot}
