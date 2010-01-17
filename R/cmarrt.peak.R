`cmarrt.peak` <-
function(cmarrt.ma,alpha,method,minrun){

peak.tmp_all=declare.bound(alpha,method,cmarrt.ma$pv.cmarrt,cmarrt.ma$pv.indep)
cmarrt.bound=indep.bound=list(Chr=NULL,Start=NULL,Stop=NULL,n.probe=NULL,min.pv=NULL,ave.pv=NULL)

if(sum(peak.tmp_all$bdd.cmarrt)>0){

for(i in 1:length(unique(cmarrt.ma$data.sort$regID))){
cregID=which(cmarrt.ma$data.sort$regID==unique(cmarrt.ma$data.sort$regID)[i])
peak.tmp=list(bdd.cmarrt=peak.tmp_all$bdd.cmarrt[cregID],bdd.indep=peak.tmp_all$bdd.indep[cregID])
if(sum(peak.tmp$bdd.cmarrt)>0){
cmarrt.tmp=rm.small.peak(peak.tmp$bdd.cmarrt,minrun,cmarrt.ma$data.sort$chr[cregID],cmarrt.ma$data.sort$start[cregID],cmarrt.ma$data.sort$stop[cregID],cmarrt.ma$pv.cmarrt[cregID])
cmarrt.bound=list(Chr=c(cmarrt.bound$Chr,cmarrt.tmp$chr.start),Start=c(cmarrt.bound$Start,cmarrt.tmp$loc.start),Stop=c(cmarrt.bound$Stop,cmarrt.tmp$loc.stop),n.probe=c(cmarrt.bound$n.probe,cmarrt.tmp$size),min.pv=c(cmarrt.bound$min.pv,cmarrt.tmp$minPv),ave.pv=c(cmarrt.bound$ave.pv,cmarrt.tmp$avePv))
}
if(sum(peak.tmp$bdd.indep)>0){
indep.tmp=rm.small.peak(peak.tmp$bdd.indep,minrun,cmarrt.ma$data.sort$chr[cregID],cmarrt.ma$data.sort$start[cregID],cmarrt.ma$data.sort$stop[cregID],cmarrt.ma$pv.indep[cregID])
indep.bound=list(Chr=c(indep.bound$Chr,indep.tmp$chr.start),Start=c(indep.bound$Start,indep.tmp$loc.start),Stop=c(indep.bound$Stop,indep.tmp$loc.stop),n.probe=c(indep.bound$n.probe,indep.tmp$size),min.pv=c(indep.bound$min.pv,indep.tmp$minPv),ave.pv=c(indep.bound$ave.pv,indep.tmp$avePv))
}
}
final=list(cmarrt.bound=cmarrt.bound,indep.bound=indep.bound)
if(length(final$cmarrt.bound$Chr)==0)print('Zero peaks found under correlation structure after postprocessing')
if(length(final$indep.bound$Chr)==0)print('Zero peaks found under independence after postprocessing')

return(final)
}

if(sum(peak.tmp_all$bdd.indep)>0&sum(peak.tmp_all$bdd.cmarrt)==0){

for(i in 1:length(unique(cmarrt.ma$data.sort$regID))){
cregID=which(cmarrt.ma$data.sort$regID==unique(cmarrt.ma$data.sort$regID)[i])
peak.tmp=list(bdd.cmarrt=peak.tmp_all$bdd.cmarrt[cregID],bdd.indep=peak.tmp_all$bdd.indep[cregID])
if(sum(peak.tmp$bdd.indep)>0){
indep.tmp=rm.small.peak(peak.tmp$bdd.indep,minrun,cmarrt.ma$data.sort$chr[cregID],cmarrt.ma$data.sort$start[cregID],cmarrt.ma$data.sort$stop[cregID],cmarrt.ma$pv.indep[cregID])
indep.bound=list(Chr=c(indep.bound$Chr,indep.tmp$chr.start),Start=c(indep.bound$Start,indep.tmp$loc.start),Stop=c(indep.bound$Stop,indep.tmp$loc.stop),n.probe=c(indep.bound$n.probe,indep.tmp$size),min.pv=c(indep.bound$min.pv,indep.tmp$minPv),ave.pv=c(indep.bound$ave.pv,indep.tmp$avePv))
}
}
final=list(indep.bound=indep.bound)
print('Zero bound probe found under correlation structure')

return(final)
}

if(sum(peak.tmp_all$bdd.indep)==0){
print('Zero bound probe found /No enrichment')
}

}

