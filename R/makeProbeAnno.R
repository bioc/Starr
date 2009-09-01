
writePosFile <- function(bpmap, file) {
	if(class(bpmap) == "character") {
		cat("Reading bpmap file: ", bpmap)	
		bpmap <- readBpmap(bpmap)
	}

	chr <- names(bpmap)
	lengths <- sapply(bpmap, function(x) {x[["seqInfo"]][["numberOfHits"]]})
	cat("Creating .pos annotation for ")
	for(i in 1:length(chr)) {
		cat(chr[i], " ")
		probe_id <- c()
		if(i == 1) {
			probe_id <- 1:lengths[chr[i]]
		}
		else {
			probe_id <- (sum(lengths[chr[1:i-1]])+1):sum(lengths[chr[1:i]])
		}
		seq_id <- rep(NA, lengths[chr[i]])
		count <- rep(NA, lengths[chr[i]])
		chromosome <- rep(chr[i], lengths[chr[i]])	
		position <- bpmap[[chr[i]]][["startpos"]]
		length <- bpmap[[chr[i]]][["probeseq"]]
		length <- sapply(length, nchar)
		
		pos <- data.frame(PROBE_ID=probe_id, SEQ_ID=seq_id, CHROMOSOME=chromosome, POSITION=position, COUNT=count, LENGTH=length, stringsAsFactors=F)
		if( i == 1 ) {
			write.table(pos, file=file, quote=FALSE, row.names = F, col.names = T, sep="\t")
		}
		else {
			write.table(pos, file=file, quote=FALSE, row.names = F, col.names = F, sep="\t", append=T)
		}
	}
	cat("\n")
}


makeProbeAnno <- function(posFile=NULL, bpmap=NULL, probeIDAsStrings=F) {
	if(is.null(bpmap) & (class(posFile) == "character")) {
		dat <- data.frame()
		if(probeIDAsStrings) {
			dat <- read.table(posFile, sep="\t", quote="", header=T, colClasses=c("character", "character", "character", "integer", "integer", "integer"))
		}
		else {
			dat <- read.table(posFile, sep="\t", quote="", header=T, colClasses=c("integer", "character", "character", "integer", "integer", "integer"))
		}
		return(posToProbeAnno(dat))
	}
	else if( (!is.null(bpmap)) & (is.null(posFile))) {
		if(class(bpmap) == "character") {
			bpmap <- readBpmap(bpmap)
		}
		chr <- names(bpmap)
		lengths <- sapply(bpmap, function(x) {x[["seqInfo"]][["numberOfHits"]]})
		nprobes <- sum(lengths)
		
		probe_id <- 1:nprobes
		if(probeIDAsStrings) {
			probe_id <- as.character(probe_id)
		}
		seq_id <- rep(NA, nprobes)
		chromosome <- unlist(sapply(chr, function(x) {rep(x, lengths[x])}))
		position <- unlist(sapply(chr, function(x) {bpmap[[x]][["startpos"]]}))
		count <- rep(NA, nprobes)
		length <- unlist(sapply(chr, function(x) {bpmap[[x]][["probeseq"]]}))
		length <- sapply(length, nchar)
		dat <- data.frame(PROBE_ID=probe_id, SEQ_ID=seq_id, CHROMOSOME=chromosome, POSITION=position, COUNT=count, LENGTH=length, stringsAsFactors=F)
		return(posToProbeAnno(dat))
	}
}




bpmapToProbeAnno <- function(bpmap, verbose=T, uniqueSeq=T) {
	if(class(bpmap) == "character") {
		cat("Reading bpmap file: ", bpmap)	
		bpmap <- readBpmap(bpmap)
	}
	
	chr <- names(bpmap)
	
	if(uniqueSeq) {
	  if(verbose) {
	    cat("Checking for uniqueness of probes on the chip\n")
	  }
	  sequences <- as.vector(unlist(lapply(bpmap, function(x) {return(x[["probeseq"]])})))
	  seq_index <- 1:length(sequences)

	  seq2index <- tapply(seq_index, INDEX=sequences, identity)
	  unique <- rep(NA, length(sequences))
	  seq2index <- lapply(seq2index, function(x) {if(length(x) == 1) {unique[x[1]] <<- 0; return(x)} else {unique[x] <<- length(x)-1; ; return(x)}})
	}

	lengths <- sapply(bpmap, function(x) {x[["seqInfo"]][["numberOfHits"]]})
	nprobes <- sum(lengths)

	probe_id <- 1:nprobes
	chromosome <- unlist(sapply(chr, function(x) {rep(x, lengths[x])}))
	probe_id <- split(probe_id, chromosome)
	position <- lapply(chr, function(x) {bpmap[[x]][["startpos"]]})
	names(position) <- chr
	count <- rep(0, nprobes)
	if(uniqueSeq) {
	  count <- unique
	}
	count <- split(count, chromosome)
	length <- lapply(chr, function(x) {nchar(bpmap[[x]][["probeseq"]]) - 1})
	names(length) <- chr


	middle_pos <- lapply(chr, function(x) {position[[x]]+length[[x]]})
	names(middle_pos) <- chr
	
	if (verbose) 
        cat("Creating probeAnno mapping for chromosome ")
	
	end <- list()
	for(i in 1:length(chr)) {
		if (verbose) 
            cat(chr[i], "")
		end[[chr[i]]] <- position[[chr[i]]]+length[[chr[i]]]
		end[[chr[i]]] <- end[[chr[i]]][order(middle_pos[[chr[i]]])]
		position[[chr[i]]] <- position[[chr[i]]][order(middle_pos[[chr[i]]])]
		probe_id[[chr[i]]] <- probe_id[[chr[i]]][order(middle_pos[[chr[i]]])]
		count[[chr[i]]] <- count[[chr[i]]][order(middle_pos[[chr[i]]])]	
	}

	thisMapping = new.env(hash = TRUE, size = length(chr) * 4)
		for (i in 1:length(chr)) {
			assign(paste(chr[i], "start", sep = "."), as.integer(position[[chr[i]]]), envir = thisMapping)
			assign(paste(chr[i], "end", sep = "."), as.integer(end[[chr[i]]]), envir = thisMapping)
			assign(paste(chr[i], "index", sep = "."), probe_id[[chr[i]]], envir = thisMapping)
			assign(paste(chr[i], "unique", sep = "."), count[[chr[i]]], envir = thisMapping)
		}
	myProbeAnno <- new("probeAnno", map = thisMapping)
	if (verbose) 
        cat("Done.\n")
	myProbeAnno
}

