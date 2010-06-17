heatmapplot <- function(profiles, colpal=c("black","dark blue","dark green", "green","gold", "yellow"), abl=NULL, subset=NULL) {

	fillIn <- function(a) {
		c <- 1
		r <- which(is.na(a))
		groups <- sapply(1:length(r), function(x) {if(x>1){if(r[x]-r[x-1] != 1) {c<<-c+1;return(c)} else{return(c)}} else{return(c)}})
		fills <- tapply(r, INDEX=groups, identity)
		lapply(fills, function(x) {if(any(x == 1)){a[x]<<-a[max(x)+1]}else if(any(x == length(a))) {a[x]<<-a[min(x)-1]} else {a[x] <<- seq(a[min(x)-1], a[max(x)+1], length=length(x)+2)[2:(length(x)+1)]}})
		return(a)
	}

	mat <- list2matrix(profiles)[["profile"]]
	ran <- range(mat, na.rm=TRUE)
	if(! is.null(subset)) {
		mat <- mat[,subset]
	}
	means <- apply(mat, 2, median, na.rm=T)
	mat <- apply(mat, 2, function(x) {if(any(is.na(x))) {fillIn(x)} else {x}})

	colfct=colorRampPalette(colpal)
	palette = colfct(30)

	levelplot(mat[, order(means)], at=seq(ran[1],ran[2], length=30), cuts = 29, col.regions = palette, aspect="fill", xlab="position", ylab="profile", panel=function(...){panel.levelplot(...);if(! is.null(abl)){panel.abline(v=abl);};panel.axis(side="bottom", at=501,labels="TSS", outside=TRUE)})
}
