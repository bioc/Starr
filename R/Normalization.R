

# wrapper for normalization methods
normalize.Probes <- function(eSet, method=NULL, ratio=FALSE, ip, control, description, fkt=median, featureData=FALSE, targets=NULL, ...) {
	
	cat(paste("Normalizing probes with method: ", method, "\n", sep=""))
	## package affy: normalize.loess
	## package limma: normalizeBetweenArrays
	## package rMAT: MAT normalization
	normalizedMatrix <- switch(method, loess = normalize.loess(exprs(eSet)), 
										none = normalizeBetweenArrays(exprs(eSet), method="none", targets=targets, ...), 
										scale = normalizeBetweenArrays(exprs(eSet), method="scale", targets=targets, ...), 
										quantile = normalizeBetweenArrays(exprs(eSet), method="quantile", targets=targets, ...), 
										Aquantile = normalizeBetweenArrays(exprs(eSet), method="Aquantile", targets=targets, ...), 
										Gquantile = normalizeBetweenArrays(exprs(eSet), method="Gquantile", targets=targets, ...), 
										Rquantile = normalizeBetweenArrays(exprs(eSet), method="Rquantile", targets=targets, ...),
										Tquantile = normalizeBetweenArrays(exprs(eSet), method="Tquantile", targets=targets, ...),
										vsn = normalizeBetweenArrays(exprs(eSet), method="vsn", targets=targets, ...),
										rankpercentile = rankPercentile.normalize(exprs(eSet)),
										substract = substract(exprs(eSet), ...))
	
	exprs(eSet) <- normalizedMatrix
	preproc(experimentData(eSet)) <- list(normalization = "method")
	if(ratio) {
		return(getRatio(eSet, ip, control, description, fkt, featureData))
	}
	else {
		return(eSet)
	}
	
}



substract <- function(matrix, fun) {
	means <- apply(matrix, 2, fun)
	matrix <- t(apply(matrix, 1, function(x) {x-means}))
	matrix
}


## rank percentile normalization
rankPercentile.normalize <- function(matrix) {

	cat("Rank Percentile normalization\n")
	nrow <- nrow(matrix)
	for(i in 1:dim(matrix)[2]) {
		matrix[,i] <- (rank(matrix[,i])-1) / (nrow-1)
	}
	
	return(matrix)
}





#expressionSet2TilingSet <- function(eSet) {
#	#cat("Converting tilingSet to ExpressionSet\n")
#	featureChromosome <- as.vector(featureData(eSet)$chr)
#	featurePosition <- featureData(eSet)$pos
#	featureCopyNumber <- as.integer(rep(1, length=length(featurePosition)))
#	exprs <- exprs(eSet)
#	featureSequence <- featureData(eSet)$seq
#	
#	
#	newSet <- new('tilingSet', featureChromosome=featureChromosome, featureSequence=featureSequence,
#				featurePosition=featurePosition, featureCopyNumber=featureCopyNumber, exprs=exprs, experimentData=experimentData(eSet))				
#}

#normalize.MAT <- function(eSet, ...) {
#	#if(verbose)
#	#	cat("Normalizing probes...\n")
#	ScSet <- expressionSet2TilingSet(eSet)
#	ScSetNorm <- NormalizeProbes(ScSet, ...)
#	mat <- exprs(ScSetNorm)
#	rownames(mat) <- featureNames(eSet)
#	colnames(mat) <- colnames(exprs(eSet))
#	mat
#}




## build ration over several arrays
getRatio <- function(eSet, ip, control, description, fkt=median, featureData=F) {
	
	ipExpressionSet <- eSet[,ip]
	controlExpressionSet <- eSet[,control]
	nrow <- dim(exprs(ipExpressionSet))[1]
	
	cat("Calculating ratio\n")
	
	exprmat <- matrix(nrow=nrow, ncol=1)
	if(dim(exprs(ipExpressionSet))[2] > 1 | dim(exprs(controlExpressionSet))[2] > 1) {
		v1 <- apply(exprs(ipExpressionSet), 1, fkt)
		v2 <- apply(exprs(controlExpressionSet), 1, fkt)
		exprmat[,1] <- v1 - v2
	}
	else {
		exprmat[,1] <- exprs(ipExpressionSet)[,1] - exprs(controlExpressionSet)[,1]
	}
	
	colnames(exprmat) <- description
	rownames(exprmat) <- rownames(exprs(eSet))

	dat2 <- data.frame(type = description)
	rownames(dat2) <- colnames(exprmat)
	pd2 <- new("AnnotatedDataFrame", data = dat2)

	cat("Building new ExpressionSet\n")
	
	if(featureData) {
		es2 <- new("ExpressionSet", exprs = exprmat, phenoData = pd2, featureData = featureData(eSet))
	}
	else {
		es2 <- new("ExpressionSet", exprs = exprmat, phenoData = pd2)
	}
	es2

}
