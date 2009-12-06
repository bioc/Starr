sortbygenomic <-
function(chr,start,stop,logR){
chr=as.character(chr)
GAP=1000 #gap in bps to be considered two segments
sort.data=matrix(0,nrow=length(chr),ncol=3)
chr.sort=rep(0,length(chr))
chrID=unique(chr)
regID=NULL
j=1
k=1

for(i in 1:length(chrID)){
sub.start=start[which(as.character(chr)==chrID[i])]
sub.stop=stop[which(as.character(chr)==chrID[i])]
reg.tmp=c(which(sort(sub.start)[-1]-sub.stop[order(sub.start)][-length(sub.stop)]>GAP),length(sub.stop))
regID=c(regID,rep(c(k:(k+length(reg.tmp)-1)),c(reg.tmp[1],(reg.tmp[-1]-reg.tmp[-length(reg.tmp)]))))
k=k+length(reg.tmp)
sub.logR=logR[which(chr==chrID[i])]
sub.chr=rep(chrID[i],length(which(chr==chrID[i])))
sort.data[j:(j+length(sub.start)-1),]=cbind(sub.start,sub.stop,sub.logR)[order(sub.start),]
chr.sort[j:(j+length(sub.start)-1)]=sub.chr
j=j+length(sub.start)
}
sort.data=data.frame(cbind(regID,chr.sort,data.frame(sort.data)))
colnames(sort.data)=c('regID','chr','start','stop','logR')
return(as.data.frame(sort.data))
}

