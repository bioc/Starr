match_ac <- function(dictionary, text, complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE) {
  #cat("Constructing search tree\n")
  result = .Call("find_ac2", dictionary=as.character(dictionary), number_of_entries=as.integer(length(dictionary)), 
			    text=as.character(text), ntext=as.integer(length(text)), complementary=as.integer(complementary),
			    reverse=as.integer(reverse), reverse_complementary=as.integer(reverse_complementary), PACKAGE="Starr")
  #cat("done.\n")
  return(result)
}

remap <- function(bpmap, nseq, bsg, chromosome_wise=FALSE, mySeqs=NULL, complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE) {

	sequences <- as.vector(unlist(lapply(bpmap, function(x) {
		    return(x[["probeseq"]])
		})))
	pmx <- as.vector(unlist(lapply(bpmap, function(x) {
		    return(x$pmx)
		})))

	pmy <- as.vector(unlist(lapply(bpmap, function(x) {
		    return(x$pmy)
		})))

	pos <- rep(NA, length(sequences))
	chr <- rep(NA, length(sequences))
	strand <- rep(NA, length(sequences))
	#occurence <- unique
	if(is.null(mySeqs)) {
		mySeqs <- seqnames(bsg)
	}



	
	steps <- seq(0,length(sequences)+nseq, by=nseq)
	steps[length(steps)] <- length(sequences)
	

	if(chromosome_wise) {
		for(i in 1:length(mySeqs)) {
			currSeq <- getSeq(bsg, strand="+", names=mySeqs[i])
				seqn <- seqnames(bsg)
			for(j in 2:length(steps)) {
				cat("Searching sequences ", steps[j-1], "to", steps[j], "in ", mySeqs[i], "\n")
		
				matches <- match_ac(sequences[(steps[j-1]+1):steps[j]], currSeq, complementary, reverse, reverse_complementary)
			myMatch <- sapply((steps[j-1]+1):steps[j], function(x) {if(matches$index[[x-steps[j-1]]] == -2) {pos[x] <<- -2; chr[x] <<- NA; strand[x] <<- NA}
										else if(matches$index[x-steps[j-1]] == -1) {}
										else {if(is.na(pos[x])) {pos[x] <<- matches$index[x-steps[j-1]]; chr[x] <<- seqn[matches$text[x-steps[j-1]]+1]; strand[x] <<- matches$strand[x-steps[j-1]]} else {pos[x] <<- -2; chr[x] <<- NA; strand[x] <<- NA}}})

		}
		}
	}
	else {
		currSeq <- getSeq(bsg, strand="+")
		seqn <- seqnames(bsg)
		for(j in 2:length(steps)) {
			#cat("Searching sequences ", steps[j-1], "to", steps[j], "\n")
		
			matches <- match_ac(sequences[(steps[j-1]+1):steps[j]], currSeq, complementary, reverse, reverse_complementary)
			myMatch <- sapply((steps[j-1]+1):steps[j], function(x) {if(matches$index[[x-steps[j-1]]] == -2) {pos[x] <<- -2; chr[x] <<- NA; strand[x] <<- NA}
										else if(matches$index[x-steps[j-1]] == -1) {}
										else {if(is.na(pos[x])) {pos[x] <<- matches$index[x-steps[j-1]]; chr[x] <<- seqn[matches$text[x-steps[j-1]]+1]; strand[x] <<- matches$strand[x-steps[j-1]]} else {pos[x] <<- -2; chr[x] <<- NA; strand[x] <<- NA}}})

		}
	}
	
	cat(round(length(which(! is.na(chr)))/length(chr)*100, digits=2), "% of the probes could be mapped uniquely.\n")

	splittedByChr <- split(1:length(chr), factor(chr))

	newbpmap <- list()
	chrnames <- names(splittedByChr)

	for(i in 1:length(splittedByChr)) {
		newbpmap[[chrnames[i]]] <- list(seqInfo=list(name=chrnames[i], groupname="", fullname=chrnames[i], version="", mapping="onlypm", number=i, numberOfHits=length(splittedByChr[[i]])), 
						pmx=pmx[splittedByChr[[i]]], pmy=pmy[splittedByChr[[i]]], mmx=NULL, mmy=NULL,
						probeseq=sequences[splittedByChr[[i]]], strand=strand[splittedByChr[[i]]], startpos=pos[splittedByChr[[i]]])
	}
	return(newbpmap)
}

