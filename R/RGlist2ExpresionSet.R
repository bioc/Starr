
RGlist2ExpresionSet <- function(RGlist) {
	ncy3 <- dim(RGlist$G)[2]
	ncy5 <- dim(RGlist$R)[2]
	mat <- cbind(RGlist$G, RGlist$R)
	rownames(mat) <- RGlist$genes$PROBE_ID
	colnames(mat) <- c(paste(colnames(mat)[1:ncy3], RGlist$targets$Cy3, sep="_"), 
						paste(colnames(mat)[(ncy3+1):(ncy3+ncy5)], RGlist$targets$Cy5, sep="_"))
	
	dat <- data.frame(file = c(RGlist$targets$FileNameCy3, RGlist$targets$FileNameCy5), 
								type = factor(c(RGlist$targets$Cy3, RGlist$targets$Cy5)))
	rownames(dat) <- colnames(mat)
	pd <- new("AnnotatedDataFrame", data = dat)
	
	
	es <- new("ExpressionSet", exprs = mat, phenoData = pd)

	es
}

