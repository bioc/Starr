
# author: Benedikt Zacher

# This function reads one or more cel files and stores it as an ExpressionSet.


readCelFile <- function(bpmap, cel_files, names, type, experimentData=NULL, featureData=T, log.it=T) {
	
	if(class(bpmap) == "character") {
		cat(paste("Reading bpmap file: ", bpmap, "\n", sep=""))
		bpmap <- readBpmap(bpmap)
	}
	
	cat("Reading cel file(s)\n")
	cel_header <- readCelHeader(cel_files[1])
	
	# Only reads PM probes, MM probes are not used for analysis
	cel_indices <- lapply(bpmap, function(x) xy2indices(x[["pmx"]], x[["pmy"]], nr = cel_header[["cols"]]))
	names(bpmap) <- names(cel_indices)

	data <- readCelIntensities(cel_files, unlist(cel_indices))
	colnames(data) <- names
	#rownames(data) <- 1:dim(data)[1]

	if(log.it) {
		cat("Taking log2 transformation of data\n")
		data <- log2(data)
	}
	
	cat("Creating phenoData\n")
	dat <- data.frame(type = factor(type), CEL = cel_files) 
	metadata <- data.frame(labelDescription = c("Description of experiment", "CEL files"), row.names = c("type", "CEL"))
	rownames(dat) <- names
	pd <- new("AnnotatedDataFrame", data = dat, varMetadata = metadata)

	## experimentData
	if(is.null(experimentData)) {
		experimentData <- new("MIAME")
	}
	if(log.it) {
		preproc(experimentData) <- list(transformation = "log", normalization = "none")
	}
	else {
		preproc(experimentData) <- list(transformation = "none", normalization = "none")
	}
	
	
	if(featureData) {
		cat("Creating featureData\n")
		fd <- new("AnnotatedDataFrame", 
			data = data.frame(
				chr = factor(unlist(lapply(as.vector(unlist(lapply(bpmap, function(x) {x[["seqInfo"]][["name"]]}))), function(y) {rep(y, length(cel_indices[[y]]))}))), 
				seq = unlist(lapply(bpmap, function(x) {x[["probeseq"]]})), 
				pos = unlist(lapply(bpmap, function(x) {x[["startpos"]]})),
				stringsAsFactors = FALSE), 
			varMetadata = data.frame(
				labelDescription = rep(c("Chromosome", "Probe sequence", "Probe start"))))
			
		cat("Creating ExpressionSet\n")
		es <- new("ExpressionSet", exprs = data, phenoData = pd, experimentData = experimentData, featureData = fd)
	}
	else {
		cat("Creating ExpressionSet\n")
		es <- new("ExpressionSet", exprs = data, phenoData = pd, experimentData = experimentData)
	}
	
	featureNames(es) <- 1:dim(data)[1]
	es
	
}

    
