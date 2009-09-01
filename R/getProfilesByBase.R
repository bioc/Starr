
getProfilesByBase <- function(eSet, probeAnno, chr, gffAnno, upstream, downstream) {
	
	
	profile <- list()
	intensity <- as.vector(exprs(eSet))
	
	cat("chromosome ")
	for(i in 1:length(chr)) {
		cat(chr[i]," ")
		index <- c(probeAnno[paste(chr[i], "index", sep=".")], probeAnno[paste(chr[i], "index", sep=".")])
		bounds <- c(probeAnno[paste(chr[i], "start", sep=".")], probeAnno[paste(chr[i], "end", sep=".")]+1) ## hier +1
		index <- index[order(bounds)]
		bounds <- sort(bounds)
		pos <- 1:max(probeAnno[paste(chr[i], "end", sep=".")]+1)
		ints <- whichIn(pos, bounds)
		
		# mapping for probe to positions in vector index
		probe2interval <- tapply(1:length(index), INDEX=index, identity)
		probe2interval <- lapply(probe2interval, function(x) {x[1]:(x[2]-1)})
		indices <- probeAnno[paste(chr[i], "index", sep=".")]

		intervals_not_covered_by_probes <- unique(ints[which(ints %in% unlist(probe2interval) == F)])
		probe2interval[["NA"]] <- intervals_not_covered_by_probes[-which(intervals_not_covered_by_probes == Inf | intervals_not_covered_by_probes == -Inf)]
		for_split <- sapply(1:length(probe2interval), function(x) {if(length(probe2interval) == x) {rep(NA, length(probe2interval[[x]]))} else{rep(indices[x], length(probe2interval[[x]]))}})
		interval2probe <- split(unlist(for_split), unlist(probe2interval))
	
		temp_int <- ints[-which(ints == Inf | ints == -Inf)]
		pos2indices <- lapply(pos, function(x) {if(ints[x] == Inf | ints[x] == -Inf) {return(NaN)} else {return(interval2probe[[ints[x]]])}})#[-which(ints == Inf | ints == -Inf)]

		
		geneAnno <- gffAnno[which(gffAnno$chr == chr[i]),]

		max_end <- max(probeAnno[paste(chr[i], "end", sep=".")]+1)
		for(j in 1:dim(geneAnno)[1]) {
			mypos <- list()
			from <- NULL
			to <- NULL
			if(geneAnno[j,]$strand == 1) {
				from <- geneAnno[j,]$start-upstream
				to <- geneAnno[j,]$end+downstream
			}
			else {
				from <- geneAnno[j,]$end+upstream
				to <- geneAnno[j,]$start-downstream
			}
			
			
			if(from < 1 & to <= max_end) {
				mypos <- pos2indices[1:to]
			}
			else if(from >= 1 & to > max_end) {
				mypos <- pos2indices[from:max_end]
			}
			else if(from < 1 & to > max_end) {
				mypos <- pos2indices[1:max_end]
			}
			else if(from >= 1 & to <= max_end){
				mypos <- pos2indices[from:to]
			}
			profile[[geneAnno[j,]$name]] <- sapply(mypos, function(x) {mean(intensity[x], na.rm=T)})
			if(from < 1) {
				profile[[geneAnno[j,]$name]] <- append(profile[[geneAnno[j,]$name]], rep(NA, abs(from)+1), 0)
			}
			if(to > max_end) {
				profile[[geneAnno[j,]$name]] <- append(profile[[geneAnno[j,]$name]], rep(NA, abs(to-max_end)+1))
			}
			if(geneAnno[j,]$start == -1) {
				profile[[geneAnno[j,]$name]] <- rev(profile[[geneAnno[j,]$name]])
			}
			
			region <- geneAnno[j,]$end - geneAnno[j,]$start
			up <- vector(mode="numeric")
			down <- vector(mode="numeric")
			
			if(upstream != 0) {
				up <- profile[[geneAnno[j,]$name]][1:upstream]
			}
			if(downstream != 0) {
				down <- profile[[geneAnno[j,]$name]][(upstream+1+region+1):(upstream+1+region+downstream)]
			}
			
			if(region == 0) {
				profile[[geneAnno[j,]$name]] <- list(upstream=up, region=profile[[geneAnno[j,]$name]][length(up)+1], downstream=down)
			}
			else {
				profile[[geneAnno[j,]$name]] <- list(upstream=up, region=profile[[geneAnno[j,]$name]][(upstream+1):(upstream+1+region)], downstream=down)
			}
			
		}
		
	}
	cat("\n")
	return(profile)
}
