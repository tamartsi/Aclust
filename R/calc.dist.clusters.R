calc.dist.clusters <-
function(clust.1, clust.2, type = "single", dist.type = "spearman"){
	## caculates the distance between two sets of vectors
	## the distance is either the min, max or mean of distances between the each of the vectors in one set
	## and each of the vectors in the other set
	## the distance function is dist.type
	stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman", "pearson", "euclid")))
	
	clust.1 <- as.matrix(clust.1)
	size.1 <- ncol(clust.1)
	if (size.1 == 1) return(calc.dist.clust.point(clust.2, clust.1, type = type, dist.type = dist.type))
	
		distances <- apply(clust.1, 2, function(x){
			calc.dist.clust.point(clust.2, x, type = type, dist = dist.type)
		})
		
	if (type == "single") return(min(distances)) else
		if (type == "complete") return(max(distances)) else
			return(mean(distances))
	
}
