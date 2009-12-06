`declare.bound` <-
function(alpha,method,pv.cmarrt,pv.indep){
pv.adj.cmarrt=p.adjust(pv.cmarrt,method)
pv.adj.indep=p.adjust(pv.indep,method)
bdd.cmarrt=rep(0,length(pv.cmarrt))
bdd.cmarrt[which(pv.adj.cmarrt<=alpha)]=1
bdd.indep=rep(0,length(pv.indep))
bdd.indep[which(pv.adj.indep<=alpha)]=1
final=list(bdd.cmarrt=bdd.cmarrt,bdd.indep=bdd.indep)
return(final)
}

