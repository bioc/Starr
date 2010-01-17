backgd.sd <-
function(regID,logR,M=NULL,chr,start,stop){
maxWin=50
N=length(logR)
tmp=density(logR)
leftlim=min(tmp$x[which(tmp$y>10e-04)])
center=tmp$x[which(tmp$y==max(tmp$y))]
yL=logR[which(logR<=center&logR>=leftlim)]
yR=2*center - yL
backgd=c(yL,yR)
sigma=sd(backgd)

if(length(M)==0)cutoff=quantile(backgd,prob=0.98) else cutoff=min(quantile(backgd,prob=0.98),quantile(logR,prob=1-M))

outlier=which(logR>cutoff)
if(outlier[1]!=1)outlier=c(1,outlier)
if(outlier[length(outlier)]!=N)outlier=c(outlier,N)
outlier=sort(unique(c(which(regID[-length(regID)]-regID[-1]!=0),outlier)))

outlierID=rep(0,N)
outlierID[outlier]=1
logR_rmout=logR
logR_rmout[outlier]=NA

getGAP=rm.small.peak(outlierID,minrun=1,chr,start,stop,pv=rep(0,N))
gapID=sort(which(getGAP$size>10))
cont_start=1
cont_end=N
if(length(gapID)>0){
if(getGAP$peak.start[gapID[1]]<=50)gatID=gapID[-1]
if(getGAP$peak.end[gapID[length(gapID)]]>=(N-50))gatID=gapID[-(length(gapID))]
cont_start=c(cont_start,(getGAP$peak.end[gapID]+1))
cont_end=c((getGAP$peak.start[gapID]-1),cont_end)
}
cont_reg=cont_end-cont_start
too_short=which(cont_reg<=10)
if(length(too_short)>0){
cont_start=cont_start[-too_short]
cont_end=cont_end[-too_short]
}

seg=cont_start

rho.tmp=matrix(0,nrow=length(seg),ncol=maxWin)
seg_length=rep(0,length(seg))
for(i in 1:length(seg)){
rho.tmp[i,]=acf(na.omit(logR_rmout[cont_start[i]:cont_end[i]]),plot=FALSE)$acf[2:(maxWin+1)]
seg_length[i]=length(na.omit(logR_rmout[cont_start[i]:cont_end[i]]))
rho.tmp[i, which(is.na(rho.tmp[i,]))]=0
}

id=which(seg_length<=200)
if(length(id)<length(seg)&length(id)>0){
if(length(id)<(length(seg)-1))rho=apply(rho.tmp[-id,],2,mean)else rho=c(rho.tmp[-id,])
if(abs(rho[1])<= 2/sqrt(mean(seg_length[-id]))){
rho=rep(0,maxWin)
print('Correlation structure is negligible')
}
}

if(length(id)==length(seg)){
print('Warning:Contiguous regions too short for reliable estimates of correlation structure')
rho=rep(0,maxWin)
}
if(length(id)==0){
rho=apply(rho.tmp,2,mean)
if(abs(rho[1])<= 2/sqrt(mean(seg_length))){
rho=rep(0,maxWin)
print('Correlation structure is negligible')
}
}

result=list(center=center,sigma=sigma,rho=rho) 
return(result)
}

