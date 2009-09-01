
\name{getRatio}
\alias{getRatio}
\title{ Building ratio over experiments }
\description{
This function calculates the ratio over experiments.
}
\usage{
getRatio(eSet, ip, control, description, fkt=median, featureData=F) 
}
\arguments{
\item{ eSet}{An ExpressionSet, containing the logged raw intensities}
\item{ip}{ a boolean or integer vector, that indicate, which columns in the matrix are IP experiments}
\item{control}{ a boolean or integer vector, that indicate, which columns in the matrix are CONTROL or REFERENCE experiments}
\item{description}{ description of the new data (e.g. IPvsCONTROL)}
\item{fkt}{ mean or median to calculate the averaged intensity over replicates}
\item{featureData}{ if TRUE, featureData is added to the new ExpressionSet }
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\examples{
##
# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))

# cels <- c(file.path(dataPath,"Rpb3_IP_chr1.cel"), file.path(dataPath,"wt_IP_chr1.cel"), 
# 	file.path(dataPath,"Rpb3_IP2_chr1.cel"))
# names <- c("rpb3_1", "wt_1","rpb3_2")
# type <- c("IP", "CONTROL", "IP")
# rpb3Chr1 <- readCelFile(bpmapChr1, cels, names, type, featureData=TRUE, log.it=TRUE)

# ips <- rpb3Chr1$type == "IP"
# controls <- rpb3Chr1$type == "CONTROL"

# rpb3_rankpercentile <- normalize.Probes(rpb3Chr1, method="rankpercentile")
# description <- c("Rpb3vsWT")
# rpb3_rankpercentile_ratio <- getRatio(rpb3_rankpercentile, ips, controls, description, fkt=median, featureData=FALSE)
}
\keyword{manip}