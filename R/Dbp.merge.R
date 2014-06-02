Dbp.merge <-
function(ordr.vec, thresh.dist, bp.thresh.dist, location.vec,  dist.type = "spearman"){
	stopifnot( is.element(tolower(dist.type), c("spearman", "pearson","euclid" )))
	
	### initialize variables
	le <- dim(ordr.vec)[2]
	dist.clust <- rep(0, le - 1)
	which.clust <- 1:le
#	ind.clust.min <- which.clust
	
	# first we look at the maximal number of neighboring distances such that their distance is smaller than thresh.dist
	max.d <- 1
	for (d in 1:min(bp.thresh.dist, length(location.vec))){
		if (min(diff(location.vec, lag = d)) <= bp.thresh.dist) max.d <- d
		else break	
	}
	
	## calculates distance vectors
	dist.mat.by.nr <- matrix(0, ncol = le, nrow = max.d) 
	for (nr in 1:max.d){
		dist.mat.by.nr[nr, ] <- calc.dist.d.neighbor(ordr.vec, nr, dist.type)
	}
	
	## how many neighbors count for each site?
	loc.diff.mat <- matrix(Inf, nrow = max.d, ncol = le)
	for (i in 1:max.d){
		diff.temp <- diff(location.vec, lag = i)
		loc.diff.mat[i, 1:length(diff.temp)] <- diff.temp
	}
	
	nr.number <- apply(loc.diff.mat, 2, function(x) max(0, which(x < bp.thresh.dist)))

		
	###  merging sites until no more to merge
	
	for (k in 1:(le - 1)){  # for each site 
		
		if (le - k > nr.number[k]) d.cur <- nr.number[k] else ## set the number of neighbors to consider to current site
			d.cur <- le - k
		
		if (d.cur > 0) {
		
		for (nr in d.cur:1){ ## start from furthest possible neighbor
			if (which.clust[k] == which.clust[k + nr]) break 
			
			cur.dist <- dist.mat.by.nr[nr,k] 
			if (cur.dist < thresh.dist ) {
				which.clust <- update.clust.indicator(which.clust, k, k + nr )
			break ## since all the closer neighbors where merged
			}
		}}  ## end if d.cur > 0
		
		
	}
	
	return(which.clust)
}
