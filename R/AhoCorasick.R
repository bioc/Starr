match_ac <- function(dictionary, text, complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE, nseq) {
  #cat("Constructing search tree\n")
  result = .Call("find_ac2", dictionary=as.character(dictionary), number_of_entries=as.integer(length(dictionary)), 
			    text=as.character(text), ntext=as.integer(length(text)), complementary=as.integer(complementary),
			    reverse=as.integer(reverse), reverse_complementary=as.integer(reverse_complementary), nseq=as.integer(nseq), PACKAGE="Starr")
  #cat("done.\n")
  return(result)
}

remap <- function(bpmap=NULL, seqs=NULL, nseq=NULL, path="", complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE, return_bpmap=FALSE) {
	
	sequences <- vector(mode="character")
	if(! is.null(bpmap)) {
	sequences <- as.vector(unlist(lapply(bpmap, function(x) {
		    return(x[["probeseq"]])
		})))
	}
	else{sequences <- seqs}
	if(is.null(nseq)) {
	  nseq=length(sequences)
	}
	#sequences <- c("ATG", "ATC", "CAT", "TG")
	
	myFASTA <- dir(path)[grep(".fa", dir(path))]
	matches <- match_ac(sequences, file.path(path, myFASTA), complementary, reverse, reverse_complementary, nseq)

	
	cat(round(length(which(matches$index >= 0))/length(matches$index)*100, digits=2), "% of the probes could be mapped uniquely.\n")
	
	if(return_bpmap) {

	  pmx <- as.vector(unlist(lapply(bpmap, function(x) {
		    return(x$pmx)
		})))

	  pmy <- as.vector(unlist(lapply(bpmap, function(x) {
		    return(x$pmy)
		})))
	  splittedByChr <- split(1:length(matches$text), factor(matches$text))

	  newbpmap <- list()
	  chrnames <- names(splittedByChr)

	  for(i in 1:length(splittedByChr)) {
		  newbpmap[[chrnames[i]]] <- list(seqInfo=list(name=chrnames[i], groupname="", fullname=chrnames[i], version="", mapping="onlypm", number=i, numberOfHits=length(splittedByChr[[i]])), 
						  pmx=pmx[splittedByChr[[i]]], pmy=pmy[splittedByChr[[i]]], mmx=NULL, mmy=NULL,
						  probeseq=sequences[splittedByChr[[i]]], strand=matches$strand[splittedByChr[[i]]], startpos=matches$index[splittedByChr[[i]]])
	  }

	  return(newbpmap)
	}
        else {
	  return(matches)
	}
}
