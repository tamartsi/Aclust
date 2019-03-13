###called by Aclust. No need to export

calc.dist.d.neighbor <-
function(ordr.vec, d, dist.type = "spearman"){
	## calucate the distances between vectors, and their d neighbor. 
	
	
	### initialize variables
	le <- dim(ordr.vec)[2]
	dist.vec.d <- rep(0, le - d)
	
	### calculate the distances vector
	if (is.element(dist.type , c("pearson", "spearman"))){
	for (i in 1:(le-d)){
		dist.vec.d[i] <- 1- abs(cor(ordr.vec[,i], ordr.vec[,i+d], use = "complete.obs", method = dist.type))
	} } else if (dist.type == "euclid"){
			inds.rm <- union(which(is.na(ordr.vec[,i])), which(is.na(ordr.vec[,i + d])))
			if (length(inds.rm) > 0) dist.vec.d[i] <- sqrt(sum((ordr.vec[-inds.rm, i] - ordr.vec[-inds.rm, i + d])^2)) else
				dist.vec.d[i] <- sqrt(sum((ordr.vec[,i] - ordr.vec[,i + d])^2))
			}
			
	dist.vec.d <- c(dist.vec.d, rep(Inf, d))		
	
	return(dist.vec.d)
}
