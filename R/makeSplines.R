
makeSplines <- function(profiles, df=1000) {
	
	stats <- c()
	l1 <- round(mean(sapply(profiles[["profile"]], function(x) {length(x[["upstream"]])})))
	l2 <- round(mean(sapply(profiles[["profile"]], function(x) {length(x[["region"]])})))
	l3 <- round(mean(sapply(profiles[["profile"]], function(x) {length(x[["downstream"]])})))
	#print(stats)
	x_pos <- function(gene) {
		bf <- length(gene[["upstream"]])
		ss <- length(gene[["region"]])
		as <- length(gene[["downstream"]])
		#print(bf)
		#print(ss)
		#print(as)
		positions <- c( seq(1, bf, length=l1), 
						seq(bf+1, bf+ss, length=l2), 
						seq(bf+ss+1, bf+ss+as, length=l3) )
		#print("here")
		positions
	}
	
	spline <- function(gene) {
					spl <- c()
					pos <- c()
					pred <- c()
					if(length(which(is.na(unlist(gene)))) >= 1) {
						pos <- 1:length(unlist(gene))
						pos <- pos[-which(is.na(unlist(gene)))]
						pred <- unlist(gene)[-which(is.na(unlist(gene)))]
					}
					else {
						pos <- 1:length(unlist(gene))
						pred <- unlist(gene)
					}
					
					#print(pos)
					#print(pred)
					
					if(length(which(is.na(unlist(gene))))/length(unlist(gene)) < 0.5) {
					#print(length(pos))
					#print(length(pred))
						spl <- predict(sm.spline(pos, pred, df=df), x_pos(gene))
					}
					else {
					#	spl <- rep(NA, length(x_pos(gene)))
					}
					
					list(upstream=spl[1:l1], region=spl[(l1+1):(l1+l2)], 
						downstream=spl[(l1+l2+1):(l1+l2+l3)])
	}

	profiles[["profile"]] <- lapply(profiles[["profile"]], spline)
	profiles
}

