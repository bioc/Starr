
# function for writing a gff file
writeGFF <- function(expressionSet, probeAnno, file) {

	chr <- unique(unlist(lapply(ls(probeAnno), strsplit, split="(\\.start|\\.unique|\\.end|\\.index)", perl=T)))
	
	source <- "Starr"
	feature <- colnames(exprs(expressionSet))
	strand <- "."
	frame <- "."
	cat("Creating gff annotation for ")
	for(i in 1:length(chr)) {
		chrom <- chr[i]
		cat("\n", chrom)
		chr_start <- probeAnno[paste(chr[i], "start", sep=".")]
		chr_end <- probeAnno[paste(chr[i], "end", sep=".")]
		chr_index <- probeAnno[paste(chr[i], "index", sep=".")]
		if(length(feature) == 1) {
			gff <- data.frame(Chromosome=chrom, Source=source, Feature=feature, Start=chr_start, Stop=chr_end, 
								ChIPsignal=exprs(expressionSet)[chr_index,], Strand=strand, Frame=frame, Attributes=chr_index)
			#gff <- gff[-(which(is.na(gff$ChIPsignal))),]
			write.table(gff, file=file, quote=FALSE, row.names = F, col.names = F, sep="\t", append=TRUE)
		}
		else {
			for(j in 1:length(feature)) {
				gff <- data.frame(Chromosome=chrom, Source=source, Feature=feature[j], Start=chr_start, Stop=chr_end, 
								ChIPsignal=exprs(expressionSet)[chr_index,j], Strand=strand, Frame=frame, Attributes=chr_index)
				#gff <- gff[-(which(is.na(gff$ChIPsignal))),]
				write.table(gff, file=file, quote=FALSE, row.names = F, col.names = F, sep="\t", append=TRUE)
			}
		}
	}
	cat("\nDone.\n")	
}
