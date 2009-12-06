cmarrt.ma <-
function(eSet, probeAnno, chr=NULL, M=NULL,frag.length,window.opt='fixed.probe'){
if(length(M)==0)M.bd=0.02 else M.bd=M

user.opt=window.opt
if(is.null(chr)) {
	chr <- unique(sapply(strsplit(ls(probeAnno), "\\."), function(x) {x[1]}))
}

start=c()
stop=c()
index=c()

start <- as.vector(sapply(chr, function(x) {probeAnno[paste(x, "start", sep=".")]}))
stop <- as.vector(sapply(chr, function(x) {probeAnno[paste(x, "end", sep=".")]}))
index <- as.vector(sapply(chr, function(x) {probeAnno[paste(x, "index", sep=".")]}))
chrom <- sapply(chr, function(x) {rep(chr, length(probeAnno[paste(x, "index", sep=".")]))})

logR=as.vector(exprs(eSet)[index,])
data=sortbygenomic(chrom,start,stop,logR)
param=backgd.sd(data$regID,data$logR,M=M.bd,data$chr,data$start,data$stop)
out=ma.stat(data$regID,data$chr,data$start,data$stop,data$logR,frag.length,
param$center,param$sigma,param$rho,window.opt=user.opt)
if(length(out)==1){
print('Terminated')
}else{
final=list(data.sort=data,ma=out$ma,z.cmarrt=out$z.cmarrt,z.indep=out$z.indep,
pv.cmarrt=out$pv.cmarrt,pv.indep=out$pv.indep)
return(final)
}
}

