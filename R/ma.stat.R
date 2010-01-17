`ma.stat` <-
function(regID,chr,start,stop,logR,frag.length,center,sigma,rho,window.opt='fixed.probe'){
    
     rho=c(rho,rep(0,100)) 
     coeff = rep(NA, length(rho))
     cmarrt.sd = rep(NA, length(rho))
     for (k in 2:(length(rho) + 1)) {
        coeff[1:(k - 1)] = 2 * (k - c(1:(k - 1)))/k
        cmarrt.sd[k - 1] = sqrt((1 + sum(coeff[1:(k - 1)] * rho[1:(k - 
            1)]))) * sigma
    }
cmarrt.sd=c(sigma,cmarrt.sd)

space=start[-1]-stop[-length(stop)]
space=space[which(space>0&space<100)]
overlap=stop[-length(stop)]-start[-1]
overlap=overlap[which(overlap>0&overlap<100)]
Length=stop-start+1
Length=Length[which(Length>0&Length<500)]
if(length(overlap)>=length(space))w.ave=round((frag.length-mean(overlap))/(mean(Length)-mean(overlap)))
if(length(overlap)<length(space))w.ave=round((frag.length+mean(space))/(mean(Length)+mean(space)))

if(w.ave>1){

uniq.reg=c(which(regID[-length(regID)]-regID[-1]!=0))+1
uniq.reg=c(1,uniq.reg,(length(regID)+1))
ma=NULL
w=NULL
ma2 = NULL
w2 = NULL
if(window.opt=='fixed.gen.dist'){
sub.mid.all = .5 * (start+stop)

ret <- .Call(
"r__entry_pt_fgd2", 
as.real(ma),
as.integer(w),
as.integer(uniq.reg),
as.integer(length(uniq.reg)),
as.real(sub.mid.all), 
as.integer(start),
as.integer(stop),
as.real(logR),
as.integer(frag.length), PACKAGE="Starr");

ma = unlist(ret[1])
w = unlist(ret[2])

}


if(window.opt=='fixed.probe'){
for(i in 1:(length(uniq.reg)-1)){
sub.start=start[uniq.reg[i]:(uniq.reg[(i+1)]-1)]
sub.stop=stop[uniq.reg[i]:(uniq.reg[(i+1)]-1)]
sub.logR=logR[uniq.reg[i]:(uniq.reg[(i+1)]-1)]
if(length(sub.logR)>w.ave){
ma.tmp=filter(sub.logR,rep(1/w.ave,w.ave),sides=2)
 for(j in 1:floor(w.ave/2))ma.tmp[j]=mean(sub.logR[1:j])
 for(j in (length(sub.logR)-floor(w.ave/2)+1):length(sub.logR))ma.tmp[j]=mean(sub.logR[j:length(ma.tmp)])
}else ma.tmp=rep(mean(sub.logR),length(sub.logR)) 
ma=c(ma,ma.tmp)
}
w=rep(w.ave,length(start))
}

uniq.w=sort(unique(w))
z.cmarrt=ma
z.indep=ma
for(i in 1:length(uniq.w)){
z.cmarrt[which(w==uniq.w[i])]=(ma[which(w==uniq.w[i])]-center)/(cmarrt.sd[uniq.w[i]]/sqrt(uniq.w[i]))
z.indep[which(w==uniq.w[i])]=(ma[which(w==uniq.w[i])]-center)/(sigma/sqrt(uniq.w[i]))
}

pv.cmarrt=pnorm(z.cmarrt,lower.tail=FALSE)
pv.indep=pnorm(z.indep,lower.tail=FALSE)

final=list(ma=ma,z.cmarrt=z.cmarrt,z.indep=z.indep,pv.cmarrt=pv.cmarrt,pv.indep=pv.indep)
return(final)
}else print('Warning:Window size=0. Cannot compute moving average. Check fragment size and array resolution--Terminated')

}

