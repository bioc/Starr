
## The following functions are used to extract the profiles (method="middle") of annotated features ##

mapFeatures <- function(probeAnno, gffAnno, upstream, downstream, chr) {
	intervals <- list()
	
	#cat("Mapping chromosome")
	for(i in 1:length(chr)) {
		intervals[[i]] <- list()
		#cat(" ", i)
		rows <- which(gffAnno$chr == chr[i])
		watson <- which(gffAnno$chr == chr[i] & gffAnno$strand == 1)
		crick <- which(gffAnno$chr == chr[i] & gffAnno$strand == -1)
		bounds <- c(gffAnno[rows,]$start, gffAnno[watson,]$start-upstream, gffAnno[crick,]$end+upstream,
					gffAnno[rows,]$end, gffAnno[watson,]$end+downstream, gffAnno[crick,]$start-downstream)
		features <- as.vector(gffAnno[rows,]$name)
		features_for_intervals <- c(as.vector(gffAnno[rows,]$name), as.vector(gffAnno[watson,]$name), as.vector(gffAnno[crick,]$name), 
					  as.vector(gffAnno[rows,]$name), as.vector(gffAnno[watson,]$name), as.vector(gffAnno[crick,]$name))
		intervals_of_features <- features_for_intervals[order(bounds)]
		features2interval <- lapply(features, function(x) {w <- which(intervals_of_features == x); 
															out <- list()
															out[["upstream"]] <- w[1]:(w[2]-1)
															out[["region"]] <- w[2]:(w[3]-1)
															out[["downstream"]] <- w[3]:(w[4]-1)
															out})
		names(features2interval) <- features
		bounds <- sort(bounds)
		length <- probeAnno[paste(chr[i], "end", sep=".")] - probeAnno[paste(chr[i], "start", sep=".")] + 1
		pos <- probeAnno[paste(chr[i], "start", sep=".")]+round(length/2)
		ints <- whichIn(pos,bounds)
		intervals[[i]][["features2interval"]] <- features2interval
		intervals[[i]][["interval2index"]] <- tapply(probeAnno[paste(chr[i], "index", sep=".")],INDEX=ints,identity)
		intervals[[i]][["interval2position"]] <- tapply(pos,INDEX=ints,identity)
	}
	#cat("\n")
	intervals
}


getFeature <- function(feature, k, mapping, expr_mat, gffAnno, reverse=T) {
	out <- list();
	
	intervals <- mapping[[k]][["features2interval"]][[feature]][["upstream"]]
	out[["upstream"]] <- as.vector(unlist(sapply(intervals, function(x) {expr_mat[mapping[[k]][["interval2index"]][[as.character(x)]],]})))
	
	
	
	intervals <- mapping[[k]][["features2interval"]][[feature]][["region"]]
	out[["region"]] <- as.vector(unlist(sapply(intervals, function(x) {expr_mat[mapping[[k]][["interval2index"]][[as.character(x)]],]})))
	
	
	intervals <- mapping[[k]][["features2interval"]][[feature]][["downstream"]]
	out[["downstream"]] <- as.vector(unlist(sapply(intervals, function(x) {expr_mat[mapping[[k]][["interval2index"]][[as.character(x)]],]})))
	
	
	if((gffAnno[which(gffAnno$name == feature),]$strand == -1)){# & reverse) {
		out[["region"]] <- rev(out[["region"]])
		temp_upstream <- out[["upstream"]]
		out[["upstream"]] <- rev(out[["downstream"]])
		out[["downstream"]] <- rev(temp_upstream)
	}
	
	out
}

getIntensities <- function(eSet, chr, mapping, gffAnno, upstream, downstream) {
		expr_mat <- exprs(eSet)
		distr <- list()
		names <- c()
		#cat("Getting probe intensities from chromosome")
		for(i in 1:length(mapping)) {
			#cat(" ", i)
			distr <- append(distr, lapply(names(mapping[[i]][["features2interval"]]), getFeature, i, mapping, expr_mat, gffAnno))
			if(length(names) == 0) {
				names <- names(mapping[[i]][["features2interval"]])
			}
			else {
				names <- append(names,names(mapping[[i]][["features2interval"]]))
			}
		}
		#cat("\n")
		names(distr) <- names
		distr
}




fill <- function(pos, val, distance, spacing, strand) {
	i <- 1
	reverse <- F
	if(strand == -1) {
		val <- rev(val)
		reverse <- T
	}
	
	while(i <= (length(pos)-1)) {
		if((pos[i]+(distance)) < pos[i+1]) {
			pos <- append(pos, pos[i]+spacing, i)
			
			val <- append(val, NA, i-1)
		}
		i <- i+1
	}
	
	if(reverse) {
		return(rev(val))
	}
	else {
		return(val)
	}
}

fillNA <- function(distribution, mapping, upstream, downstream, gffAnno, distance, spacing) {

	distr <- distribution

	for(j in 1:length(mapping)) {
		features <- names(mapping[[j]][["features2interval"]])
		take_from_gff <- which(gffAnno$name %in% features)
		for( i in 1:length(features)) {
			
			if(gffAnno[take_from_gff[i],]$strand == 1) {
				
				if(upstream > 0) {
					intervals <- mapping[[j]][["features2interval"]][[features[i]]][["upstream"]]
					positions <- as.vector(unlist(sapply(intervals, function(x) {mapping[[j]][["interval2position"]][[as.character(x)]]})))
					values <- distr[[features[i]]][["upstream"]]
					
					distr[[features[i]]][["upstream"]] <- fill(c(gffAnno[take_from_gff[i],]$start-upstream,positions, gffAnno[take_from_gff[i],]$start), values, distance, spacing, 1)
				}
				else {
					distr[[features[i]]][["upstream"]] <- vector(mode="numeric")
				}
			
				intervals <- mapping[[j]][["features2interval"]][[features[i]]][["region"]]
				positions <- as.vector(unlist(sapply(intervals, function(x) {mapping[[j]][["interval2position"]][[as.character(x)]]})))
				values <- distr[[features[i]]][["region"]]
				
				distr[[features[i]]][["region"]] <- fill(c(gffAnno[take_from_gff[i],]$start,positions, gffAnno[take_from_gff[i],]$end), values, distance, spacing, 1)
				
				if(downstream > 0) {
					intervals <- mapping[[j]][["features2interval"]][[features[i]]][["downstream"]]
					positions <- as.vector(unlist(sapply(intervals, function(x) {mapping[[j]][["interval2position"]][[as.character(x)]]})))
					values <- distr[[features[i]]][["downstream"]]
					
					distr[[features[i]]][["downstream"]] <- fill(c(gffAnno[take_from_gff[i],]$end,positions, gffAnno[take_from_gff[i],]$end+downstream), values, distance, spacing, 1)
				}
				else {
					distr[[features[i]]][["downstream"]] <- vector(mode="numeric")
				}
			}
			else if(gffAnno[take_from_gff[i],]$strand == -1) {
				if(upstream > 0) {
					intervals <- mapping[[j]][["features2interval"]][[features[i]]][["downstream"]]
					positions <- as.vector(unlist(sapply(intervals, function(x) {mapping[[j]][["interval2position"]][[as.character(x)]]})))
					values <- distr[[features[i]]][["upstream"]]
					
					distr[[features[i]]][["upstream"]] <- fill(c(gffAnno[take_from_gff[i],]$end,positions, gffAnno[take_from_gff[i],]$end+upstream), values, distance, spacing, -1)
				}
				else {
					distr[[features[i]]][["upstream"]] <- vector(mode="numeric")
				}
				
				intervals <- mapping[[j]][["features2interval"]][[features[i]]][["region"]]
				positions <- as.vector(unlist(sapply(intervals, function(x) {mapping[[j]][["interval2position"]][[as.character(x)]]})))
				values <- distr[[features[i]]][["region"]]
				
				distr[[features[i]]][["region"]] <- fill(c(gffAnno[take_from_gff[i],]$end,positions, gffAnno[take_from_gff[i],]$start), values, distance, spacing, -1)
				
				if(downstream > 0) {
					intervals <- mapping[[j]][["features2interval"]][[features[i]]][["upstream"]]
					positions <- as.vector(unlist(sapply(intervals, function(x) {mapping[[j]][["interval2position"]][[as.character(x)]]})))
					values <- distr[[features[i]]][["downstream"]]
					
					distr[[features[i]]][["downstream"]] <- fill(c(gffAnno[take_from_gff[i],]$start-downstream,positions, gffAnno[take_from_gff[i],]$start), values, distance, spacing, -1)
				}
				else {
					distr[[features[i]]][["downstream"]] <- vector(mode="numeric")
				}
			}
		}
	}
	
	distr
}

#prof <- getProfiles(rpb3_rankpercentile_ratio, probeAnno, geneAnno[1:117,], 0, 0, feature="ORF", borderNames=c("start", "stop"), method="middle", sameLength=F)
#plot(unlist(prof[["profile"]][["YAR075W"]]), ylim=c(-0.3,0.7))
#abline(h=0)


sameLength <- function(distribution, method=c("upstream", "downstream")) {
	
	upstream <- min(sapply(distribution, function(x) {length(x[["upstream"]])}))
	downstream <- min(sapply(distribution, function(x) {length(x[["downstream"]])}))
	for(i in 1:length(distribution)) {
		if("upstream" %in% method) {
			distribution[[i]][["upstream"]] <- distribution[[i]][["upstream"]][(length(distribution[[i]][["upstream"]])-upstream+1):length(distribution[[i]][["upstream"]])]
		}
		if("downstream" %in% method) {
			distribution[[i]][["downstream"]] <- distribution[[i]][["downstream"]][1:downstream]
		}
	}
	
	distribution
}
