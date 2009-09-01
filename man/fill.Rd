
\name{fill}
\alias{fill}
\title{Fill large spaces in one profile with NA}
\description{
This function fills spaces in a profiles, that exceed a given distance.
It is called by fillNA.
}
\usage{
fill(pos, val, distance, spacing, strand)
}
\arguments{
\item{pos}{a vector with genomic position of the signals}
\item{val}{a vector, with ChIP signals}
\item{distance}{spaces, exceeding this distance are filled with NA}
\item{spacing}{probe spacing on the chip}
\item{strand}{watson (+1), or crick (-1)}
}
\value{
a vector of signals, filled with NA
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{ internal }
