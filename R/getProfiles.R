getProfiles <- function(eSet, probeAnno, gffAnno, upstream, downstream, feature="ORF", borderNames, method, sameLength=T, fill=T, distance=8, spacing=4) {
	
	
	chr <- unique(unlist(lapply(ls(probeAnno), strsplit, split="(\\.start|\\.unique|\\.end|\\.index)", perl=T)))
	chr <- intersect(chr, gffAnno$chr)
	
	profile <- list()
	profile[["ID"]] <- colnames(exprs(eSet))
	profile[["upstream"]] <- upstream
	profile[["downstream"]] <- downstream
	profile[["method"]] <- method
	profile[["borderNames"]] <- borderNames
	profile[["feature"]] <- feature
	
	if(method == "middle") {
		cat("Mapping annotated features to spotted probes \n")
		featureMap <- mapFeatures(probeAnno, gffAnno, upstream, downstream, chr)
		cat("Getting probe intensities \n")
		distribution <- getIntensities(eSet, chr, featureMap, gffAnno)
		#distribution <- lapply(distribution, function(x) {if(! all(is.na(distribution))) {return(x)}})
		if(fill) {
			cat("Filling gaps with NAs\n")
			distr <- fillNA(distribution, featureMap, upstream, downstream, gffAnno, distance, spacing)
			
			if(sameLength) {
				cat("Making equal lengths for upstream/downstream regions\n")
				distr <- sameLength(distr)
				profile[["profile"]] <- distr
				return(profile)
			}
			else {
				profile[["profile"]] <- distr
				return(profile)
			}
		}
		else {
			profile[["profile"]] <- distribution
			return(profile)
		}
	}
	else if(method == "basewise") {
		profile[["profile"]] <- getProfilesByBase(eSet, probeAnno, chr, gffAnno, upstream, downstream)
		return(profile)
	}
}