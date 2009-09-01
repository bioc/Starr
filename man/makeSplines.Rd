
\name{makeSplines}
\alias{makeSplines}
\title{Fit splines to profiles}
\description{
This function uses the pspline package to fit spilnes to each entry in a list of profiles.
}
\usage{
makeSplines(profiles, df=1000)
}
\arguments{
\item{profiles}{a list as it is created by the getProfiles package.}
\item{df}{the degree of freedom of the fit}
}
\value{
a list as it is created by the getProfiles function.
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[pspline]{sm.spline}}, \code{\link[pspline]{predict}}}

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

# probeAnnoChr1 <- bpmapToProbeAnno(bpmapChr1)
# transcriptAnno <- read.gffAnno(file.path(dataPath, "transcriptAnno.gff"), feature="transcript")

# profile <- getProfiles(rpb3_rankpercentile_ratio, probeAnnoChr1, transcriptAnno, 500, 500, feature="transcript", borderNames=c("TSS", "TTS"), method="basewise", sameLength=T, fill=T, distance=8, spacing=4)
# profile_splines <- makeSplines(profile)
}
\keyword{manip}
