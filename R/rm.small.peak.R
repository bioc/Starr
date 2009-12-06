`rm.small.peak` <-
function(bdd.method,minrun,chr,start,stop,pv){

peak.start=NULL
peak.end=NULL
size=NULL

if(length(bdd.method)>minrun&length(which(bdd.method==0))>0){
ret <- .Call("r__rm_small_peak", as.integer(bdd.method), length(bdd.method), as.double(minrun), PACKAGE="Starr")

peak.start = unlist(ret[1])
peak.end = unlist(ret[2])
size = unlist(ret[3])
}

if(length(which(bdd.method==1))==length(bdd.method)&length(bdd.method)>minrun){
peak.start=1
peak.end=length(bdd.method)
size=peak.end-peak.start+1
}

if(length(peak.start)>0){
chr.start=as.character(chr[peak.start])
chr.stop=as.character(chr[peak.end])
loc.start=start[peak.start]
loc.stop=stop[peak.end]
avePv=rep(0,length(loc.stop))
minPv=rep(0,length(loc.stop))

for(i in 1:length(loc.stop)){
avePv[i]=mean(pv[peak.start[i]:peak.end[i]])
minPv[i]=min(pv[peak.start[i]:peak.end[i]])
}
}else{

peak.start = NULL
peak.end = NULL
size = NULL
chr.start=NULL
chr.stop=NULL
loc.start=NULL
loc.stop=NULL
avePv=NULL
minPv=NULL

}

final=list(bdd.method=bdd.method,size=size,peak.start=peak.start,peak.end=peak.end,chr.start=chr.start,chr.stop=chr.stop,loc.start=loc.start,loc.stop=loc.stop,minPv=minPv,avePv=avePv)
return(final)
}

