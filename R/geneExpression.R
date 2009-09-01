

expressionByFeature <- function(eSet, fkt, method="median") {
	features <- unlist(mget(featureNames(eSet), fkt))
	affy_ids <- names(features)
	exprByORF <- list()
	for(i in 1:length(features)) {
		nr <- which(features == features[i])
		if(method == "median") {
			exprByORF[[features[i]]] <- median(exprs(eSet)[affy_ids[nr],1], na.rm=T)
		}
		else if(method == "mean") {
			exprByORF[[features[i]]] <- mean(exprs(eSet)[affy_ids[nr],1], na.rm=T)
		}
	}
	
	unlist(exprByORF)
}


intersection <- function(means, expression) {
	names <- intersect(names(means), names(expression))
	if(length(which(is.na(names))) > 0) {
		names <- names[-which(is.na(names))]
	}
	mat <- matrix(nrow=length(names), ncol=2)
	mat[,1] <- means[names]
	mat[,2] <- expression[names]
	rownames(mat) <- names
	
	mat
}


filterGenes <- function(gffAnno, distance_us=500, distance_ds=500, minLength=-Inf, maxLength=Inf) {
	
	chr <- unique(as.vector(gffAnno$chr))
	
	order <- unlist(sapply(chr, function(x) {temp <- order(gffAnno[which(as.vector(gffAnno$chr) == x),]$start); which(as.vector(gffAnno$chr) == x)[temp]}))
	gffAnno <- gffAnno[order,]
	
	features <- c()
	chr <- "none"
	
	for(i in 1:(dim(gffAnno)[1]-1)) {
		geneLength <- as.numeric(gffAnno[i,]$end) - as.numeric(gffAnno[i,]$start)
		if(as.vector(gffAnno[i,]$chr) != chr) {
			if((as.numeric(gffAnno[i,]$end)+distance_ds<as.numeric(gffAnno[i+1,]$start)) & geneLength > minLength & geneLength < maxLength) {
				features[length(features)+1] <- gffAnno[i,]$name
			}
			
			chr <- as.vector(gffAnno[i,]$chr)
		}
		else {
			if((as.numeric(gffAnno[i,]$start)-distance_us > as.numeric(gffAnno[i-1,]$end)) & (as.numeric(gffAnno[i,]$end)+distance_ds<as.numeric(gffAnno[i+1,]$start)) & geneLength > minLength & geneLength < maxLength) {
				features[length(features)+1] <- gffAnno[i,]$name
			}
		}
	}
	
	features
}

