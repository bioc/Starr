whichIn <- function(pos,bounds) {
	a <- c(-Inf,bounds,Inf,pos)
	b <- c(1:(length(bounds)+2),rep(0,length(pos)))
	h1 <- order(a)
	hr <- b[h1]
	intpos <- cummax(hr)[which(hr==0)]-1
	intpos[intpos==0] <- -Inf
	intpos[intpos==length(bounds)] <- Inf
	intpos
}
