
list2matrix <- function(profiles) {
	
	lengths <- sapply(profiles[["profile"]], function(x) {length(unlist(x))})
	
	stopifnot(length(unique(lengths)) == 1)
	
	mat <- matrix(unlist(lapply(profiles[["profile"]], function(x) {as.vector(unlist(x))})), ncol=length(profiles[["profile"]]))
	colnames(mat) <- names(profiles[["profile"]])
	profiles[["profile"]] <- mat
	profiles
}

