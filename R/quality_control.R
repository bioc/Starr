

plotImage <- function(cel) {
	cel <- readCel(cel, readXY = TRUE)
	levelplot(log2(cel[["intensities"]])~cel[["x"]]*cel[["y"]], cuts = 16, col.regions=rainbow(22)[1:17], xlab="Y position", ylab="X position",region=T)
}


plotMA <- function(eSet, ip=NULL, control=NULL, col=NULL) {
	
	if(is.null(col)) {
		col <- "#00000040"
	}	
	
	if(inherits(eSet, "ExpressionSet")) {
		eSet <- exprs(eSet)
	}

	xy <- windowxy(length(which(ip == T)) * length(which(control == T)))
	par(mfcol=c(xy[1], xy[2]))
	names <- colnames(eSet)
	#combinations <- list()
	cat("Plotting ")
	for(i in 1:length(ip)) {
		for(j in 1:length(control)) {
			if(ip[i] == T & control[j] == T) {
				cat(paste(names[i], names[j], sep="vs"), " ")
				#combinations[[length(combinations)+1]] <- c(i, j)
				M <- eSet[,i] - eSet[,j]
				A <- (eSet[,i] + eSet[,j])/2
				
				ma.plot(A, M, pch=".", cex=0.8, col=col, main=paste(names[i], names[j], sep="vs"))
			}
		}
	}
	cat("\n")
	
}



plotBoxes <- function(eSet, col=NULL) {
	
	if(inherits(eSet, "ExpressionSet")) {
		eSet <- exprs(eSet)
	}
	
	cols <- c("#666666", "#999966", "#669966", "#996666", "#996699", "#669699", "#666999", "#999696", "#696966")
	
	if(is.null(col) & length(cols) < dim(eSet)[2]) {
		boxplot(eSet, names=NA, varwidth=T)
		mtext(side=1, at=1:dim(eSet)[2], colnames(eSet), padj=1, font=2, line=2)
	}
	else if(is.null(col) & length(cols) >= dim(eSet)[2]) {
		boxplot(eSet, col=cols, names=NA, varwidth=T)
		mtext(side=1, at=1:dim(eSet)[2], colnames(eSet), padj=1, font=2, line=2)
	}
	else if(!is.null(col)) {
		boxplot(eSet, col=col, names=NA, varwidth=T)
		mtext(side=1, at=1:dim(eSet)[2], colnames(eSet), padj=1, font=2, line=2)
	}
	
}


plotDensity <- function(eSet, oneDevice=T, main="") {
	
	
	if(inherits(eSet, "ExpressionSet")) {
		eSet <- exprs(eSet)
	}
	
	cols <- rainbow(dim(eSet)[2])

	expr <- NULL
	if(class(eSet) == "matrix") {
		expr <- eSet
	}
	else {
		expr <- eSet
	}
	densities <- list()
	ymax <- c()
	ymin <- c()
	for(i in 1:dim(eSet)[2]) {
			densities[[i]] <- density(expr[,i])
			ymax[i] <- max(densities[[i]]$y)
			ymin[i] <- min(densities[[i]]$y)
	}
	
	if(oneDevice) {
		plot(densities[[1]]$x,  densities[[1]]$y, col=cols[1], lwd=2, ylim=c(min(ymin), max(ymax)), xlim=c(min(expr), max(expr)), main=main, type="l", xlab="intensity", ylab="density")
		for(i in 2:dim(eSet)[2]) {
			lines(densities[[i]], col=cols[i], lwd=2)
		}
		legend("topleft", colnames(eSet), col=cols, lty=rep(1, dim(eSet)[2]) , lwd=rep(2, dim(eSet)[2]))
	}
	else {
		par(mfrow=windowxy(dim(eSet)[2]))
		for(i in 1:dim(eSet)[2]) {
			plot(densities[[1]], col=cols[i], lwd=2, main=colnames(expr)[i])
		}
	}
}



plotGCbias <- function (intensity, sequence, main = "") 
{
    gc <- sapply(strsplit(sequence, ""), function(x) sum(x %in% 
        c("C", "G")))
    boxplot(split(intensity-mean(intensity, na.rm=TRUE), gc), main = main, varwidth = TRUE, 
        xlab = "GC bases", ylab = "intensity")
}

plotPosBias <- function (intensity, sequence, main = "", ylim = NULL) 
{
    m <- mean(intensity, na.rm=TRUE)
    means <- list()
    bases <- c("A", "T", "C", "G")
    seq <- lapply(sequence, function(x) {
        unlist(strsplit(x, ""))
    })
    cat("Calculating mean intensities for ")
    indices <- list()
    for (i in 1:length(bases)) {
        cat(bases[i], " ")
        pos <- lapply(seq, function(x) {
            which(x == bases[i])
        })
        for_split <- sapply(1:length(pos), function(x) {
            rep(x, length(pos[[x]]))
        })
        indices <- split(unlist(for_split), unlist(pos))
        means[[bases[i]]] <- sapply(indices, function(x) {
            mean(intensity[x], na.rm = T)-m
        })
    }
    cat("\n")
    if (is.null(ylim)) {
        plot(means[["A"]], pch = "A", col = "green", ylim = c(min(unlist(means)), 
            max(unlist(means))), xlab = "position in sequence", 
            ylab = "intensity", main = main)
    }
    else {
        plot(means[["A"]], pch = "A", col = "green", ylim = ylim, 
            xlab = "position in sequence", ylab = "intensity", 
            main = main)
    }
    points(means[["T"]], pch = "T", col = "blue")
    points(means[["G"]], pch = "G", col = "orange")
    points(means[["C"]], pch = "C", col = "red")
}




plotScatter <- function(eSet, density=F, cluster=T, sample=NULL, cex=1) {
	if(inherits(eSet, "ExpressionSet")) {
		eSet <- exprs(eSet)
	}
	
	if(! is.null(sample)) {
		eSet <- eSet[sample(1:nrow(eSet), sample),]
	}
	
	if(cluster) {
		correlation <- cor(eSet, use="pairwise.complete.obs")
		dist <- dist(correlation)
		clust <- hclust(dist, method="average")
		eSet <- eSet[,clust$labels[clust$order]]
	}
	
	if(density) {
		pairs(eSet, lower.panel=function(x,y) {text(mean(sum(range(x, na.rm=T)))/2, mean(sum(range(y, na.rm=T)))/2, labels=round(cor(x,y, use="pairwise.complete.obs"), 2), cex=1.5)}, upper.panel=function(x,y) {densityscatter(x,y,add=T,cex=cex)})
	}
	else {
		pairs(eSet, lower.panel=function(x,y) {text(mean(sum(range(x, na.rm=T)))/2, mean(sum(range(y, na.rm=T)))/2, labels=round(cor(x,y, use="pairwise.complete.obs"), 2), cex=1.5)}, upper.panel=function(x, y) {points(x,y,pch=".");  lines(lowess(x, y), col="red")})#, pch=".", lwd=1)#, col="#00000040")
	}
}



plotRatioScatter <- function(eSet, ip, control, density=F, sample=NULL, cluster=T, cex=1) {
	
	if(inherits(eSet, "ExpressionSet")) {
		eSet <- exprs(eSet)
	}
	mat <- matrix(nrow=nrow(eSet), ncol=length(which(ip == T)) * length(which(control == T)))
	names <- colnames(eSet)
	counter <- 1
	colnames <- c()
	for(i in 1:length(ip)) {
		for(j in 1:length(control)) {
			if(ip[i] == T & control[j] == T) {
				colnames[counter] <- paste(names[i], names[j], sep="vs")
				mat[,counter] <- eSet[,i]-eSet[,j]
				counter <- counter + 1
			}
		}
	}
	colnames(mat) <- colnames
	if(! is.null(sample)) {
		mat <- mat[sample(1:nrow(mat), sample),]
	}
	
	if(cluster) {
		correlation <- cor(mat, use="pairwise.complete.obs")
		dist <- dist(correlation)
		clust <- hclust(dist, method="average")
		mat <- mat[,clust$labels[clust$order]]
	}
	
	if(density) {
		pairs(mat, lower.panel=function(x,y) {text(mean(sum(range(x, na.rm=T)))/2, mean(sum(range(y, na.rm=T)))/2, labels=round(cor(x,y, use="pairwise.complete.obs"), 2), cex=1.5)}, upper.panel=function(x,y) {densityscatter(x,y,add=T,cex=cex)})
	}
	else{
		pairs(mat, lower.panel=function(x,y) {text(mean(sum(range(x, na.rm=T)))/2, mean(sum(range(y, na.rm=T)))/2, labels=round(cor(x,y, use="pairwise.complete.obs"), 2), cex=1.5)}, upper.panel=function(x, y) {points(x,y,pch=".");  lines(lowess(x, y), col="red")})
	}
}



kde2dplot <- function(x,y,                # a 2d density computed by kde2D
                      grid=50, ncol=30,nlevels=10,main=""    # the number of colors to use
                  ){
d = kde2d(x,y,n=grid)
z   <- d$z
nrz <- nrow(z)
ncz <- ncol(z)

couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)

par(mar=c(2,2,2,2))
image(d,col=couleurs,main=main)
contour(d,add=T,nlevels=nlevels)
box()
}


densityscatter = function(x,y,pch=19,cex=1,ncol=30,grid=100,
			palette="heat", add=F,...){

getfrommat = function(a){d$z[a[1],a[2]]}
palette = switch(palette,
	heat = c("grey","dark blue","red","orange","gold"),
	red = c("#E0C6C6","red"),
	green = c("#C6E0C6","green"),
	blue = c("#C6C6E0","blue"),
	black = c("#E5E5E5","black"),
	mountain = c("light green","dark green","black","dark grey","#F0F0F0"),
	palette)
todiscrete = function(t,tmin,tmax,bins){
	erg = round((t-tmin)/(tmax-tmin)*bins+0.5)
	erg = pmin(pmax(erg,1),bins)
	return(erg)
	}

colfct=colorRampPalette(palette)
colpal = colfct(ncol)
d <- kde2d(x,y,n=grid)

xdiscrete = todiscrete(x,min(x),max(x),bins=grid)
ydiscrete = todiscrete(y,min(y),max(y),bins=grid)
heatvec = unlist(apply(cbind(xdiscrete,ydiscrete),1,getfrommat))
coldiscrete = todiscrete(heatvec,min(d$z),max(d$z),bins=ncol)
if(add) {
	points(x,y,col=colpal[coldiscrete],pch=pch,cex=cex,...)
}
else {
	plot(x,y,col=colpal[coldiscrete],pch=pch,cex=cex,...)
}
}

