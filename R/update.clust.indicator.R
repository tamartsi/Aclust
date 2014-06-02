update.clust.indicator <-
function(which.clust, ind.1, ind.2){
	## merge all points between ind.1 and ind.2, and previous/later point if belong to the same cluster. 
	## which.clust is a vector of cluster numbers, so that for each point it says which cluster in belongs to. 
	## ind.1 and ind.2 are two points that should be merged into the same vector (as are all points between them)
	
	stopifnot(ind.2 > ind.1)
	
	ind.2 <- last(which(which.clust == which.clust[ind.2]))
	which.clust[(ind.1 + 1):ind.2] <- which.clust[ind.1]
	return(which.clust)
}
