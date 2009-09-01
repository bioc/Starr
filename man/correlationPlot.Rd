
\name{correlationPlot}
\alias{correlationPlot}
\title{correlation of ChIP signals to other data}
\description{
\code{correlationPlot} The correlationPlot can be used to vizualize e.g. the correlation between the mean binding intensity of specific
 regions around annotated features and gene expression. The regions around the annotated features, that should be analyzed are defined
in a data frame. Each row represents one region. In the upper panel of the plot, the correlation is plotted in a barplot. In the lower
panel, the annotated feature and the regions defined in the data frame are shown.
}
\usage{
correlationPlot(regions, labels=c("start", "stop"), ...)
}
\arguments{
\item{ regions }{ a data frame, containing four columns. Every row defines one region to be analyzed and is plotted in the lower panel. pos=start, upstream=500 and downstream=500 mean characterize the region 500 bp upstream and downstream around the start of the feature. The pos columns is a character with values out of c("start", "region", "end"). upstream and downstream ar integers, indicating how many bases upstream and downstream from the specified position in the feature are included. level is an integer, that says at which level the rectangle in the lower device should be plotted. The numeration goes from the bottom to the ceiling. cor is the correlation of the region, which is plotted in the upper panel.   }
\item{ labels }{ a character vector which holds the names of the borders of the annotated region. (e.g. c("TSS", "TTS") for transcripts)  }
\item{...}{parameters, that are passed to barplot (for plotting the upper panel)}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[graphics]{barplot}}}
\examples{
## Constructing an example data frame
pos <- c("start", "start", "start", "region", "region","region","region", "stop","stop","stop")
upstream <- c(500, 0, 250, 0, 0, 500, 500, 500, 0, 250)
downstream <- c(0, 500, 250, 0, 500, 0, 500, 0, 500, 250)
level <- c(1, 1, 2, 3, 4, 5, 6, 1, 1, 2)
cor <- seq(-1,1, length=10)
info <- data.frame(pos=pos, upstream=upstream, downstream=downstream, level=level, cor=cor, stringsAsFactors=FALSE) 
rownames(info) <- letters[1:10]

## Generate plot
correlationPlot(info)
}
\keyword{hplot}
