
\name{profileplot}
\alias{profileplot}
\title{Vizualize clusters}
\description{
Visualization of a set of ``profiles'' (i.e. a consecutive series of measurements like a time series, or the DNA binding levels along different positions on a gene). The profiles are given as the rows of a (samples x positions) matrix that contains the measurements.
Instead of plotting a line for each profile (row of the matrix), the q-quantiles for each position (column of the matrix) are calculated, where q runs through a set of representative quantiles. 
Then for each q, a line of q-quantiles is plotted along the positions. Color coding of the quantile profiles aids the interpretation of the plot: There is a color gradient from the median profile to the 0 (=min) resp. 1(=max) quantile. 
}
\usage{
profileplot(cluster, label=NULL,	at=NULL, main = "", xlim=NULL, xlab = "", xaxt = "s",xlabels = NULL, las = 1, ylim=NULL, ylab = "", fromto = c(0.05,0.95),colpal = "heat", 	nrcolors = 25, outer.col="light grey", add.quartiles = T, add = F, separate = T)
}
\arguments{
\item{cluster}{a (samples x columns) matrix with numerical entries. Each sample row is understood as a consecutive series of measurements. Missing values are not allowed so far}
\item{label}{if multiple clusters should be plotted in one diagram, the cluster labels for each item are given in this vector}
\item{at}{optional vector of length ncol(cluster), default = 1:ncol(cluster). Specifies the x-values at which the positions will be plotted.}
\item{main}{the title of the plot, standard graphics parameter}
\item{xlim}{xlimits, standard graphics parameter}
\item{xlab}{x-axis legend, standard graphics parameter}
\item{xaxt}{should an x axis be plotted at all? (="n" if not), standard graphics parameter}
\item{xlabels}{character vector. If specified, this text will be added at the ``at``-positions as x-axis labels.}
\item{las}{direction of the xlabels text. las=1: horizontal text, las=2: vertical text}
\item{ylim}{ylimits, standard graphics parameter}
\item{ylab}{y-axis legend, standard graphics parameter}
\item{fromto}{determines the smallest and the largest quantile that are plotted in colors, more distant values are plotted as outliers}
\item{colpal}{either "red","green","blue" (predefined standard color palettes in profileplot), or a vector of colors to be used instead.}
\item{nrcolors}{not very important. How many colors will the color palette contain? Usually, the default = 25 is sufficient}
\item{outer.col}{color of the outlier lines, default = "light grey". For no outliers, choose outer.col="none"}
\item{add.quartiles}{should the quartile lines be plotted (grey/black)? default=TRUE}
\item{add}{should the profile plot be added to the current plot? Defaults to FALSE}
\item{separate}{should each cluster, be plotted in a separate window? Defaults to TRUE}
}
\examples{
sampls = 100
probes = 63
at = (-31:31)*14
clus = matrix(rnorm(probes*sampls,sd=1),ncol=probes)
clus= rbind( t(t(clus)+sin(1:probes/10))+1:nrow(clus)/sampls , t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/sampls )
labs = paste("cluster",kmeans(clus,4)$cluster)

profileplot(clus,main="All data",fromto=c(0,1))
profileplot(clus,label=labs,main="Clustered data",colpal=c("heat","blue","red","topo"),add.quartiles=FALSE)
profileplot(clus,main="Same data, 4 clusters in one plot\n color gradient fromto = c(0.4,0.6), no outliers plotted",label=labs,separate=FALSE,xaxt="n",fromto=c(0.4,0.6),
	colpal=c("heat","blue","red","green"),outer.col="none")
}
\author{ Achim Tresch, Benedikt Zacher \email{tresch@lmb.uni-muenchen.de}}
\keyword{hplot}
