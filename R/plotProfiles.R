plotProfiles <- function(profiles, mfcol=NULL, mfrow=NULL, ylab="intensity", xlab="position", histograms=NULL, cluster, profileplot=T, meanprofile=T, ...) {
	
	ncluster <- length(unique(cluster))
	prof_mat <- list2matrix(profiles)
	
	## Determining partitioning of window
	if(!is.null(mfcol) | !is.null(mfrow)) {
		if(! is.null(mfcol)) {
			par(mfcol=mfcol)
		}
		else {
			par(mfrow=mfrow)
		}
	}
	else {
		myPlots <- length(which(c(profileplot, meanprofile)) == T)*ncluster
		if(! is.null(histograms)) {
			myPlots <- myPlots+length(histograms)
		}
		par <- windowxy(myPlots)
		par(mfcol=par)

	}
	cols <- c(heat.colors(1), rainbow(1,start=0.7,end=0.1) ,colorRampPalette(c("red","yellow","green"))(1), colorRampPalette(c("dark blue","blue","#7390EE","light blue"))(1),
			"red" =  colorRampPalette(c( "#940000","#A50000","#FF5C5C","#FFB9B9"))(1),"green" = colorRampPalette(c("dark green","#009700","green","#C0F5D0"))(1))

			
	## Determining xlabels and its positions ##
	xlabels <- c(paste("-", profiles[["upstream"]], sep=""), profiles[["borderNames"]], paste("+", profiles[["downstream"]], sep=""))
	at <- c()
	if(length(xlabels) == 4) {
		l1 <- sapply(profiles[["profile"]], function(x) {length(x[["upstream"]])})
		l2 <- sapply(profiles[["profile"]], function(x) {length(c(x[["upstream"]], x[["region"]]))})
		l3 <- sapply(profiles[["profile"]], function(x) {length(unlist(x))})
		at <-  c(0, mean(l1, na.rm=T), mean(l2, na.rm=T), mean(l3, na.rm=T))
	}
	else if(length(xlabels) == 3) {
		at <- c(0, mean(sapply(profiles[["profile"]], function(x) {length(x[["upstream"]])}), na.rm=T), mean(sapply(profiles[["profile"]], function(x) {length(unlist(x))}), na.rm=T))
	}
	
	
	## profileplot ##
	labs <- cluster[colnames(prof_mat[["profile"]])]
	labs <- paste("cluster", labs)
	#x <- apply(prof_mat[["profile"]], 2, function(x) {x <- x-mean(x); x <- x/sqrt(sum(x^2)); return(x)})
	#if(profileplot) {
	#	profileplot(t(prof_mat[["profile"]]),label=labs,colpal=c("heat","blue", "topo", "jamaica", "red", "green"),add.quartiles=T, xaxt="n", xlab=xlab, ylab=ylab)
	#}
	cols <- rainbow(length(unique(cluster)))
	
	## means of clusters ##
	if(meanprofile) {
	clustersets = split(1:ncol(prof_mat[["profile"]]),factor(labs))
	means <- list()
	for(i in 1:length(clustersets)) {
		mat <- prof_mat[["profile"]][,clustersets[[i]]]
		means[[i]] <- apply(mat, 1, mean, na.rm=T)
	}
	for(i in 1:length(clustersets)) {
		plot(means[[i]], xaxt="n", main=paste("cluster ", i, ", size=", length(clustersets[[i]]), sep=""), col=cols[i], ylab=ylab, xlab=xlab, ...)		
		if(length(at) == 3) {
			abline(v=at[2], lty=3)
		}
		else if(length(at) == 4) {
			abline(v=at[2], lty=3)
			abline(v=at[3], lty=3)
		}
		axis(1, labels = xlabels, at=at)
		profileplot(t(prof_mat[["profile"]][,clustersets[[i]]]),label=rep(paste("cluster", i), length(clustersets[[i]])),colpal=c("heat","blue", "topo", "jamaica", "red", "green"),add.quartiles=T, xaxt="n", xlab=xlab, ylab=ylab)
		if(length(at) == 3) {
			abline(v=at[2], lty=3)
		}
		else if(length(at) == 4) {
			abline(v=at[2], lty=3)
			abline(v=at[3], lty=3)
		}		
		axis(1, labels = xlabels, at=at)	
	}
	}


	## 	  Density	 ##
	if(! is.null(histograms)) {
	densities <- list()
	ymax <- list()
	ymin <- list()
	xmax <- list()
	xmin <- list()
	for(i in 1:length(histograms)) {
		densities[[names(histograms)[i]]] <- list()
		ymax[[i]] <- vector(mode="numeric")
		ymin[[i]] <- vector(mode="numeric")
		xmax[[i]] <- vector(mode="numeric")
		xmin[[i]] <- vector(mode="numeric")
		for(j in 1:ncluster) {
			cluster_j <- histograms[[i]][names(which(cluster == j))]
			densities[[names(histograms)[i]]][[j]] <- density(cluster_j, na.rm=T)
			ymax[[i]][j] <- max(densities[[i]][[j]]$y, na.rm=T)
			ymin[[i]][j] <- min(densities[[i]][[j]]$y)
			xmax[[i]][j] <- max(densities[[i]][[j]]$x)
			xmin[[i]][j] <- min(densities[[i]][[j]]$x)
		}
	}

	

		## Adding rugs to density plots
		#rugpos <- seq(1.5, ncluster/2+1, by=0.5)+1
		for(i in 1:length(densities)) {
			for(j in 1:length(densities[[i]])) {
				if(j == 1) {
					plot( densities[[i]][[j]], xlim=c( c(min(xmin[[i]], na.rm=T), max(xmax[[i]], na.rm=T) )), ylim=c( c(min(ymin[[i]], na.rm=T), max(ymax[[i]], na.rm=T) )), xlab="" , lwd=2, main=names(histograms)[i], col=cols[j])
		#			rug(histograms[[i]][names(which(cluster == j))], col=cols[j], ticksize=0.01, line=rugpos[j])
				}
				else {
					lines( densities[[i]][[j]], col=cols[j], lwd=2)
		#			rug(histograms[[i]][names(which(cluster == j))], col=cols[j], ticksize=0.01, line=rugpos[j])
				}
			}
			legend("topright", paste("cluster", 1:ncluster), lwd=rep(2, ncluster), col=cols[1:ncluster], bty="n")
		}
		
		#}
	}
}

