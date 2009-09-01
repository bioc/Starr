# cluster: a #lines x #at matrix containing intensity values to be plotted
# at: a (sorted) vector of length #at containing the at of the respective probes in the cluster matrix

singleclusterplot = function(
	cluster,  # an items x columns matrix with numerical entries. each item will define a line in the clusterplot
	label=NULL,	# if multiple clusters should be plotted in one diagram, the cluster labels for each item are given in this vector
	at=NULL, # at which x-positions will the columns occur in the plot?
	main = "",  # the title of the plot
	xlim=NULL,  # xlimits, standard graphics parameter
	xlab = "",  # x-axis legend, standard graphics parameter 
	xaxt = "s", # should an x axis be plotted at all? (="n" if not)
	xlabels = NULL, las = 1, # text added as x-axis labels (las=1: horizontal text, la2=2: vertical text)
	ylim=NULL,  # ylimits, standard graphics parameter
	ylab = "",	# y-axis legend, standard graphics parameter
	fromto = c(0.05,0.95),
	colpal = "heat",  # either "red","green","blue" (standard colors), or
						# a vector of colors that can be used instead of a standard color palette
	nrcolors = 25, # how many colors will the color palette contain?
	outer.col="light grey", # color of the outlier lines
	add.quartiles = T # should the quartile lines be plotted (grey/black)?
	)
{
	title(main)
	drawline = function(y,col="black",lwd=1,lty=1)
		{lines(at[1:length(y)],y,type="l",col=col,lwd=lwd,lty=lty)}
	if (outer.col!="none") apply(cluster,1,drawline,col=outer.col)
	probes = ncol(cluster)
	
	firstquart = ceiling(nrcolors/2)
	nrcolpal = 2*firstquart+2
	nrquants = 4*firstquart+3
	if (length(colpal)==1){
	colpal = switch(colpal,
		"heat" = heat.colors(nrcolors),
		"topo" = rainbow(nrcolors,start=0.7,end=0.1) ,
		"jamaica" = colorRampPalette(c("red","yellow","green"))(nrcolors),
		"blue" = colorRampPalette(c("dark blue","blue","#7390EE","light blue"))(nrcolors),
		"red" =  colorRampPalette(c( "#940000","#A50000","#FF5C5C","#FFB9B9"))(nrcolors),
		"green" = colorRampPalette(c("dark green","#009700","green","#C0F5D0"))(nrcolors)
		)
	} else 	 colpal = colorRampPalette(colpal)(nrcolors)
	colpal = colpal[c(nrcolpal:1,2:nrcolpal)]

		qline = apply(cluster,2,quantile, na.rm=T,probs=seq(fromto[1],fromto[2],length=nrquants))
	for (j in 1:(nrquants-1))
	{
	polygon(at[c(1:probes,probes:1)],
		c(qline[j,],qline[j+1,probes:1]),
		col = colpal[j],lty=0)
	}
	if (add.quartiles){
		drawline(qline[2*firstquart+2,],col="black",lwd=1)
		drawline(qline[firstquart+1,],col="grey",lwd=1)
		drawline(qline[3*firstquart+3,],col="grey",lwd=1)
	}
}

windowxy <- function(windows=1)
{
 a = floor(sqrt(windows))
 if (a^2==windows) return(c(a,a))
 if ( windows <= a*(a+1) ) return(c(a,a+1))
 return(c(a,a+2))
}



profileplot = function(
	cluster,  # an items x columns matrix with numerical entries. each item will define a line in the clusterplot
	label=NULL,	# if multiple clusters should be plotted in one diagram, the cluster labels for each item are given in this vector
	at=NULL, # at which x-positions will the columns occur in the plot?
	main = "",  # the title of the plot
	xlim=NULL,  # xlimits, standard graphics parameter
	xlab = "",  # x-axis legend, standard graphics parameter 
	xaxt = "s", # should an x axis be plotted at all? (="n" if not)
	xlabels = NULL, las = 1, # text added as x-axis labels (las=1: horizontal text, la2=2: vertical text)
	ylim=NULL,  # ylimits, standard graphics parameter
	ylab = "",	# y-axis legend, standard graphics parameter
	fromto = c(0.05,0.95),
	colpal = "heat",  # either "red","green","blue" (standard colors), or
						# a vector of colors that can be used instead of a standard color palette
	nrcolors = 25, # how many colors will the color palette contain?
	outer.col="light grey", # color of the outlier lines
	add.quartiles = T, # should the quartile lines be plotted (grey/black)?
	add = F, # should plot be added to a current plot?
	separate = T # should clusters be plotted separately?
	)
{
	if (is.null(at)) at=1:ncol(cluster)
	probes = length(at)
	if (is.null(xlim)) xlim=c(min(at),max(at))
	maxp=xlim[2]; minp=xlim[1]
	if (is.null(ylim)) ylim=c(min(cluster, na.rm=T),max(cluster, na.rm=T))

	if (is.null(label)){	
		if(add==F){
			#x11()
			plot(xlim,ylim,
			type="n",
			xaxt="n",
			xlab=xlab,
			ylab=ylab,
			main = main)
			if (xaxt!="n") axis(side=1,las=las,at=at,labels=xlabels)
		}
		singleclusterplot(cluster=cluster,label=label,at=at,main=main,xlim=xlim,xlab=xlab,xaxt=xaxt,
					xlabels=xlabels,ylim=ylim,ylab=ylab,fromto=fromto,colpal=colpal,
					nrcolors=nrcolors,outer.col=outer.col,add.quartiles=add.quartiles)

	} else
	{
	clusternames = sort(unique(label))
	nrclusters = length(clusternames)
	clustersets = split(1:nrow(cluster),factor(label))
	if (!is.list(colpal)) colpal = as.list(colpal)
	if (length(colpal)<nrclusters) colpal = rep(colpal,nrclusters)	
	if (length(main)>0) main = paste(main,"\n",sep="")
	
	if (separate==T) add=F 
	if (add==FALSE){
		#x11()
		#plot(xlim,ylim,
		#	type="n",
		#	xaxt="n",
		#	xlab=xlab,
		#	ylab=ylab)
		#if (xaxt!="n") axis(side=1,las=las,at=at,labels=xlabels)
	}
	#if (separate==T) par(mfrow=windowxy(nrclusters)) 

	for (j in seq(clusternames)){
		if (separate==T){
			plot(xlim,ylim,
				type="n",
				xaxt="n",
				xlab=xlab,
				ylab=ylab)
			if (xaxt!="n") axis(side=1,las=las,at=at,labels=xlabels)
		}
		titel = if(separate==F) main else paste(main,clusternames[j],sep="")
		singleclusterplot(cluster=cluster[clustersets[[j]],,drop=F],label=NULL,at=at,
							main=titel,
							xlim=xlim,xlab=xlab,xaxt=xaxt,
							xlabels=xlabels,ylim=ylim,ylab=ylab,fromto=fromto,colpal=colpal[[j]],
							nrcolors=nrcolors,outer.col=outer.col,add.quartiles=add.quartiles)
	} # end for
	} # end if label
}





