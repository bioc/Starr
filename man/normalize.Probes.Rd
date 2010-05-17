
\name{normalize.Probes}
\alias{normalize.Probes}
\title{ Normalization of probes}
\description{
Normalization of probe intensities with a given method.
}
\usage{
normalize.Probes(eSet, method=NULL, ratio=FALSE, ip, control, description, fkt=median, featureData=FALSE, targets=NULL, arrays = NULL, ...)
}
\arguments{
\item{eSet}{An ExpressionSet, containing the logged raw intensities}
\item{method}{character string specifying the normalization method to be used. Choices are 
				"none", "scale", "quantile", "Aquantile", "Gquantile", "Rquantile",
				"Tquantile", "vsn", "rankpercentile", "loess", "substract".}
\item{ratio}{if TRUE, the ratios are calcualted.}
\item{ip}{a boolean vector, indicating which sample are IP experiments}
\item{control}{a boolean vector, indicating which sample are CONTROL experiments}
\item{description}{description of the normalized data}
\item{fkt}{function to chose for averaging over replicates}
\item{featureData}{should the featureData of eSet be passed to the new ExpressionSet?}
\item{targets}{vector, factor or matrix of length twice the number of arrays, used to indicate target groups if method="Tquantile"}
\item{arrays}{Subset of experiments (colnames in ExpressionSet) in the ExpressionSet, that are supposed to be normalized seperately.}
\item{...}{arguments, that should be passed to the normalization methods.}
}
\details{
 The procedure calls different functions from this package or from affy and limma, depending on the method.
 \describe{
   \item{none}{Calls \code{normalizeWithinArrays} with \code{method="none"}  from package limma.}   
	 \item{scale}{Calls \code{normalizeWithinArrays} with
     \code{method="scale"}  from package limma.}   
	 \item{quantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="quantile"}  from package limma.}
   \item{Gquantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="Gquantile"}  from package limma.}
   \item{Rquantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="Rquantile"}  from package limma.}
	\item{Tquantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="Tquantile"}  from package limma.}
	 \item{Rquantile}{Calls \code{normalizeBetweenArrays} with
     \code{method="Rquantile"}  from package limma.}
	 \item{vsn}{Calls \code{normalizeBetweenArrays} with
     \code{method="vsn"}  from package limma.}
	 \item{loess}{Calls \code{normalize.loess} from package affy.}
   \item{rankpercentile}{Calls \code{rankPercentile.normalize} from this package.}
    \item{substract}{Calls \code{substract} from this package.}

 }
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[limma:normalizebetweenarrays]{normalizeBetweenArrays}}, \code{\link[affy]{normalize.loess}}
		, \code{\link[Starr]{substract}}, \code{\link[Starr]{rankPercentile.normalize}}}
\examples{
##
# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))

# cels <- c(file.path(dataPath,"Rpb3_IP_chr1.cel"), file.path(dataPath,"wt_IP_chr1.cel"), 
# 	file.path(dataPath,"Rpb3_IP2_chr1.cel"))
# names <- c("rpb3_1", "wt_1","rpb3_2")
# type <- c("IP", "CONTROL", "IP")
# rpb3Chr1 <- readCelFile(bpmapChr1, cels, names, type, featureData=TRUE, log.it=TRUE)

# rpb3_rankpercentile <- normalize.Probes(rpb3Chr1, method="rankpercentile")
}
\keyword{manip}
