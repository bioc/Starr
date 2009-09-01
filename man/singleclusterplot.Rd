
\name{singleclusterplot}
\alias{singleclusterplot}
\title{single cluster plot}
\description{
This function is called by clusterplot
}
\usage{
singleclusterplot(cluster,label=NULL,at=NULL,main = "",xlim=NULL,xlab = "", xaxt = "s",	xlabels = NULL, las = 1, ylim=NULL, ylab = "",fromto = c(0.05,0.95),colpal = "heat", nrcolors = 25,	outer.col="light grey",	add.quartiles = T)
}
\arguments{
\item{cluster}{an items x columns matrix with numerical entries. each item will define a line in the clusterplot}
\item{label}{if multiple clusters should be plotted in one diagram, the cluster labels for each item are given in this vector}
\item{at}{at which x-positions will the columns occur in the plot?}
\item{main}{the title of the plot}
\item{xlim}{xlimits, standard graphics parameter}
\item{xlab}{x-axis legend, standard graphics parameter }
\item{xaxt}{should an x axis be plotted at all? (="n" if not)}
\item{xlabels}{text added as x-axis labels}
\item{las}{las=1: horizontal text, la2=2: vertical text}
\item{ylim}{ylimits, standard graphics parameter}
\item{ylab}{y-axis legend, standard graphics parameter}
\item{fromto}{}
\item{colpal}{either "red","green","blue" (standard colors), or a vector of colors that can be used instead of a standard color palette.}
\item{nrcolors}{how many colors will the color palette contain?}
\item{outer.col}{color of the outlier lines}
\item{add.quartiles}{should the quartile lines be plotted (grey/black)?}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
