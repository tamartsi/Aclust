calc.dist.clust.point <-
function(clust, pnt, type = "single", dist.type = "spearman"){
	## caculates the distance between a set of vectors and a single vector. 
	## the distance is either the min, max or mean of distances between the vector and each of the vectors in the set
	## the distance function is dist.type
	
	stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman", "pearson", "euclid")), dim(clust)[1] == length(pnt))
	
	clust <- as.matrix(clust)
	size <- ncol(clust)
	if (size == 1){
		if (is.element(dist.type , c("pearson", "spearman"))) return(1 - abs(cor(clust, pnt, use = "complete.obs", method = dist.type))) else
	 		if (dist.type == "euclid") {
					inds.rm <- union(which(is.na(clust)), which(is.na(pnt)))
				 if (length(inds.rm) > 0) return(sqrt(sum((clust[-inds.rm] - pnt[-inds.rm])^2))) else
					 return(sqrt(sum((clust - pnt)^2)))	
			}	
	
	} else{
	
	distances <- apply(clust, 2, function(x){
		if (is.element(dist.type , c("pearson", "spearman"))) return(1 - abs(cor(x, pnt, use = "complete.obs", method = dist.type))) else
			if (dist.type == "euclid") {
				inds.rm <- union(which(is.na(x)), which(is.na(pnt)))
				if (length(inds.rm) > 0) return(sqrt(sum((x[-inds.rm] - pnt[-inds.rm])^2))) else
					return(sqrt(sum((x - pnt)^2)))	
			}	
	} )
	
	if (type == "single") return(min(distances)) else
		if (type == "complete") return(max(distances)) else
			return(mean(distances)) }
}
