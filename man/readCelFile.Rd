
\name{readCelFile}
\alias{readCelFile}
\title{ Read raw intensities from CEL files}
\description{
Function to read the raw intensities of the perfect match probes (PM) of Affymetrix CEL files into an ExpressionSet.
This function is used to read one-color data. For two-color data use the functions from the Ringo package.
}
\usage{
readCelFile(bpmap, cel_files, names, type, experimentData=NULL, featureData=T, log.it=T)
}
\arguments{
\item{ bpmap}{ 	Either a list, created by the function readBpmap() from the affy package, or the path to the bpmap file. }
\item{ cel_files }{ 	a character vector, specifying the path to the CEL files }
\item{ names }{ 	a character vector, containing the names of the experiments }
\item{ type }{ 	a character vector, containing the type of experiment, e.g. "IP" for an Immunoprecipitation, or "CONTROL" for a control or reference experiment was done}
\item{ experimentData }{	This must be an object of type MIAME, which details information about e.g., the investigator or lab where the experiment was done, an overall title, and other notes }
\item{ featureData }{ 	If TRUE, a featureData object is added to the ExpressionSet, containing information about the chromosome, position in the genome and sequence of the features }
\item{ log.it }{ 	If TRUE, logged intesities are read}
}
\value{
Returns raw intensity values in form of an ExpressionSet with additional information:
\item{assayData }{This object contains the measured probe intensities.}
\item{phenoData }{contains further description of the experiments, such as names or type}
\item{featureData }{containing information about the chromosome, position in the genome and sequence of the features}
\item{experimentData}{details information about e.g., the investigator or lab where the experiment was done}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[affxparser]{readCelIntensities}}, \code{\link[affy]{xy2indices}}}
\examples{
##
# dataPath <- system.file("extdata", package="Starr")
# bpmapChr1 <- readBpmap(file.path(dataPath, "Scerevisiae_tlg_chr1.bpmap"))

# cels <- c(file.path(dataPath,"Rpb3_IP_chr1.cel"), file.path(dataPath,"wt_IP_chr1.cel"), 
# 	file.path(dataPath,"Rpb3_IP2_chr1.cel"))
# names <- c("rpb3_1", "wt_1","rpb3_2")
# type <- c("IP", "CONTROL", "IP")
# rpb3Chr1 <- readCelFile(bpmapChr1, cels, names, type, featureData=TRUE, log.it=TRUE)
}
\keyword{IO}