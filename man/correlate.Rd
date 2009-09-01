
\name{correlate}
\alias{correlate}
\title{Correlate the values of two named vectors}
\description{
\code{correlate} calculates the correlation of a list of named vectors and a named vector. Only entries, existing in both vectors are used.
}
\usage{
correlate(regions=NULL, means, expression, method="spearman")
}
\arguments{
\item{means}{a list of named integer vectors.}
\item{expression}{a named integer vector}
\item{method}{method for calculating the correlation}
\item{ regions }{ a data frame, containing four columns. The pos columns is a character with values out of c("start", "region", "end"). upstream and downstream ar integers, indicating how many bases upstream and downstream from the specified position in the feature are included. level is an integer, that says at which level the rectangle in the lower device should be plotted. The numeration goes from the bottom to the ceiling. cor is the correlation of the region, which is plotted in the upper panel.   }
}
\value{
a list with the correlations of each element of the list to the expression vector.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
