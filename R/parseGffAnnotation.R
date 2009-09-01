
read.gffAnno <- function(gffFile, feature=NULL) {
	gff <- read.table(gffFile, sep="\t", as.is=TRUE, quote="", stringsAsFactors=F,
	header=FALSE, colClasses=c("character", "character", "character", "integer", "integer", "character", "character", "character", "character"))
	colnames(gff) = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
	stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
	
	watson <- which(gff$strand == "+" | gff$strand == "1")
	crick <- which(gff$strand == "-" | gff$strand == "-1")
	strand <- rep(NA, dim(gff)[1])
	strand[c(watson, crick)] <- c(rep(1, length(watson)), rep(-1, length(crick)))
	gff$strand <- strand
	
	split <- strsplit(gff$attributes, split=";", fixed = TRUE)
     	gff$name <- sapply(split, function(x) {
         names_split <- strsplit(x, split = "=", fixed = TRUE)
         index_of_name = match("ID", sapply(names_split, function(y) {y[1]}))
         if (!is.na(index_of_name)) {
             out <- names_split[[index_of_name]][2]
         }
         else {
             out <- "NA"
         }
         return(out)
     	})	

	if(is.null(feature)) {
		return(gff)
	}
	else {
		take <- which(gff$feature == feature)
		stopifnot(length(take) > 0)
		return(gff[take,])
	}
}
