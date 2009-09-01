
getMeans <- function(eSet, probeAnno, geneAnno, regions) {
	watson <- which(geneAnno$strand == 1)
	crick <- which(geneAnno$strand == -1)
	#correlation <- c()
	means <- list()
	for(i in 1:dim(regions)[1]) {
		anno <- geneAnno
		if(regions[i,]$pos == "start") {
			anno[watson,]$end <- anno[watson,]$start
			anno[crick,]$start <- anno[crick,]$end
		}
		else if(regions[i,]$pos == "stop") {
			anno[watson,]$start <- anno[watson,]$end
			anno[crick,]$end <- anno[crick,]$start
		}

		d <- getProfiles(eSet, probeAnno, gffAnno=anno, upstream=regions[i,]$upstream, downstream=regions[i,]$downstream, feature="unknown", borderNames="unknown", method="middle", sameLength=F, fill=F)

		means[[i]] <- sapply(d[["profile"]], function(x) {mean(unlist(x), na.rm=T)})
		#mat <- intersection(means, exprByORF)
		
		#correlation[length(correlation)+1] <- cor(mat[,1], mat[,2], use="pairwise.complete.obs")
		#mean[length(mean)+1] <- mean(mat[,1],na.rm=T)
	}
	
	#regions$mean <- mean
	#regions$cor <- correlation
	#regions
	means
}



#pos <- c("start", "start", "start", "region", "region","region","region", "stop","stop","stop")
#upstream <- c(500, 0, 250, 0, 0, 500, 500, 500, 0, 250)
#downstream <- c(0, 500, 250, 0, 500, 0, 500, 0, 500, 250)
#info <- data.frame(pos=pos, upstream=upstream, downstream=downstream, stringsAsFactors=F)


#means_rpb3 <- getMeans(normalizedExpressionSet, probeAnno, gffAnno, info)

correlate <- function(regions=NULL, means, expression, method="spearman") {

	correlation <- sapply(means, function(x) {
									mat <- intersection(x, expression)
									cor(mat[,1], mat[,2], use="pairwise.complete.obs", method=method)
								})	
	
	if(is.null(regions)) {
		return(correlation)
	}
	else {
		regions$cor <- correlation
		return(regions)
	}
}

#info_with_cor <- correlate(info, means_rpb3, exprByORF)



correlationPlot <- function(regions, labels=c("start", "stop"), ...) {

	stopifnot(all(regions$pos %in% c("start")) | all(regions$pos %in% c("stop")) | all(c("start", "region", "stop") %in% regions$pos))

	
	cols <- colorRampPalette(brewer.pal(9, "GnBu"))(dim(regions)[1])
	names <- rownames(regions)
	names(cols) <- names[order(regions$cor)]
	
	max_us_start <- c()
	max_ds_start <- c()
	max_us_stop <- c()
	max_ds_stop <- c()
	if("start" %in% regions$pos) {
		max_us_start <- max(regions[which(regions$pos == "start"),]$upstream)
		max_ds_start <- max(regions[which(regions$pos == "start"),]$downstream)
	}
	if("stop" %in% regions$pos) {
		max_us_stop <- max(regions[which(regions$pos == "stop"),]$upstream)
		max_ds_stop <- max(regions[which(regions$pos == "stop"),]$downstream)
	}
	
	

	
	if("start" %in% regions$pos & "stop" %in% regions$pos & "region" %in% regions$pos) {
		#stopifnot(all(regions$pos %in% c("start", "region", "stop")))
		par(mfrow=c(2,1), mar=c(2,3,2,1), oma=c(2,0,0,2))
		barplot(sort(regions$cor), col=cols, names.arg=names[order(regions$cor)], ...)
		length_of_lower <- max_us_start+max_ds_stop+2*max_us_stop+2*max_ds_start+max_ds_stop
		plot(c(1,length_of_lower), c(0,0), type='l', yaxt="n", xaxt="n", xlab="position relative to start/stop codon", ylim=c(0,max(regions$level)/10), lwd=2, col="white", ylab="")
		
		stop_pos <- max_us_start+max_ds_stop+2*max_us_stop+2*max_ds_start
		for(i in 1:dim(regions)[1]) {
			if(as.vector(regions[i,]$pos) == "start") {
				left <- max_us_start-regions[i,]$upstream
				right <- max_us_start+regions[i,]$downstream
				rect(left,(regions[i,]$level-1)/10, right, regions[i,]$level/10, col=cols[names[i]])
				text((left+right)/2,(regions[i,]$level/10)-0.05, names[i])
			}
			else if(as.vector(regions[i,]$pos) == "stop") {
				left <- stop_pos-regions[i,]$upstream
				right <- stop_pos+regions[i,]$downstream
				rect(left,(regions[i,]$level-1)/10, right, regions[i,]$level/10, col=cols[names[i]])
				text((left+right)/2,(regions[i,]$level/10)-0.05, names[i])
			}
			else if(as.vector(regions[i,]$pos) == "region") {
				left <- max_us_start-regions[i,]$upstream
				right <- stop_pos+regions[i,]$downstream
				rect(left,(regions[i,]$level-1)/10, right, regions[i,]$level/10, col=cols[names[i]])
				text((left+right)/2,(regions[i,]$level/10)-0.05, names[i])
			}
		}
		axis(1, at=c(0,max_us_start,max_us_start+max_ds_stop,max_us_start+max_ds_stop+2*max_us_stop+max_ds_start,max_us_start+max_ds_stop+2*max_us_stop+2*max_ds_start,max_us_start+max_ds_stop+2*max_us_stop+2*max_ds_start+max_ds_stop), 
				labels=c(paste("-", max_us_start, sep=""), labels[1], paste("+", max_ds_start, sep=""), paste("-", max_us_stop, sep=""), labels[2], paste("+", max_ds_stop, sep="")))
		abline(v=max_us_start, lty=3)
		abline(v=stop_pos, lty=3)
	}
	else if("start" %in% regions$pos) {
		#stopifnot(all(regions$pos %in% c("start")))
		par(mfrow=c(2,1), mar=c(2,3,2,1), oma=c(2,0,0,2))
		barplot(sort(regions$cor), col=cols, names.arg=names[order(regions$cor)], ...)
		
		length_of_lower <- max_us_start+max_ds_start
		plot(c(1,length_of_lower), c(0,0), type='l', yaxt="n", xaxt="n", xlab="position relative to start/stop codon", ylim=c(0,max(regions$level)/10), lwd=2, col="white", ylab="")
		for(i in 1:dim(regions)[1]) {
			if(as.vector(regions[i,]$pos) == "start") {
				left <- max_us_start-regions[i,]$upstream
				right <- max_us_start+regions[i,]$downstream
				rect(left,(regions[i,]$level-1)/10, right, regions[i,]$level/10, col=cols[names[i]])
				text((left+right)/2,(regions[i,]$level/10)-0.05, names[i])
			}
		}
		axis(1, at=c(0,max_us_start,max_us_start+max_ds_start), labels=c(paste("-", max_us_start, sep=""), labels[1], paste("+", max_ds_start, sep="")))
		abline(v=max_us_start, lty=3)
	}
	else if("stop" %in% regions$pos & !("region" %in% regions$pos)) {
		#stopifnot(all(regions$pos %in% c("stop")))
		par(mfrow=c(2,1), mar=c(2,3,2,1), oma=c(2,0,0,2))
		barplot(sort(regions$cor), col=cols, names.arg=names[order(regions$cor)], ...)
		
		length_of_lower <- max_us_stop+max_ds_stop
		plot(c(1,length_of_lower), c(0,0), type='l', yaxt="n", xaxt="n", xlab="position relative to stop/stop codon", ylim=c(0,max(regions$level)/10), lwd=2, col="white", ylab="")
		for(i in 1:dim(regions)[1]) {
			if(as.vector(regions[i,]$pos) == "stop") {
				left <- max_us_stop-regions[i,]$upstream
				right <- max_us_stop+regions[i,]$downstream
				rect(left,(regions[i,]$level-1)/10, right, regions[i,]$level/10, col=cols[names[i]])
				text((left+right)/2,(regions[i,]$level/10)-0.05, names[i])
			}
		}
		axis(1, at=c(0,max_us_stop,max_us_stop+max_ds_stop), labels=c(paste("-", max_us_stop, sep=""), labels[1], paste("+", max_ds_stop, sep="")))
		abline(v=max_us_stop, lty=3)
	}
}

#level <- c(1, 1, 2, 3, 4, 5, 6, 1, 1, 2)
#info$level <- level
#correlationPlot(info[8:10,])
