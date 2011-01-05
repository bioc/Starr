writeWIG = function (expressionSet, probeAnno, file, chr=NULL, probeLength=NULL) 
{

    stopifnot(dim(exprs(expressionSet))[2] == 1)
    if(is.null(chr)) {
    	chr <- unique(unlist(lapply(ls(probeAnno), strsplit, split = "(\\.start|\\.unique|\\.end|\\.index)", 
        	perl = T)))
    }
    source <- "Starr"
    feature <- colnames(exprs(expressionSet))
       

    cat("Creating *.wig file for:")
    for (i in 1:length(chr)) {
	chr_start <- probeAnno[paste(chr[i], "start", sep = ".")]
        chr_end <- probeAnno[paste(chr[i], "end", sep = ".")]
        chr_index <- probeAnno[paste(chr[i], "index", sep = ".")]
	if(is.null(probeLength)) {
		probeLength = round(mean(chr_end-chr_start))
    	} 
	header = paste("track type=wiggle_0 name=\"variableStep\" description=", feature, " graphType=points\n", "variableStep\t", "chrom=", chr[i],  "\tspan=", probeLength, "", sep="")

	write.table(header, file = file, quote = FALSE, row.names = F, col.names = F, sep = "\t", append = TRUE)
        chrom <- chr[i]
        cat("\n", chrom)
       
        gff <- data.frame(Start = chr_start, ChIPsignal = exprs(expressionSet)[chr_index,])
        write.table(gff, file = file, quote = FALSE, row.names = F, col.names = F, sep = "\t", append = TRUE)
       
    }
    cat("\nDone.\n")
}

