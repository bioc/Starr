plotcmarrt <- function(cmarrt.ma, freq=FALSE){
	par(mfrow=c(2,2))
	hist(cmarrt.ma$pv.indep,main='Under independence',freq=freq,xlab='p-values')
	hist(cmarrt.ma$pv.cmarrt,main='Under correlation',freq=freq,xlab='p-values')
	qqnorm(cmarrt.ma$z.indep,main='Under independence',pch='.')
	abline(0,1)
	qqnorm(cmarrt.ma$z.cmarrt,main='Under correlation',pch='.')
	abline(0,1)
}

