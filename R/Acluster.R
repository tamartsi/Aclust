#' Title A-clustering
#' An algorithm for clustering of adjacent clusters
#' @param ordr.vec An (n by m) matrix of n measurements of m random variables (methylations). The columns are ordered. 
#' @param thresh.dist A distance threshold. Two neighboring clusters are merged to a single cluster if the distance between them is above thresh.dist.
#' @param which.clust A vector of initial clusters assignments of the m variables. If it is not provided, it is taken that each site is a unique cluster. 
#' @param location.vec An m vector specifying the locations (e.g. chromosomal locations) of the variables measured in the matrix ordr.vec. 
#' @param max.dist Optional maximum length between neighboring variables permitting to cluster them together. 
#' @param type Type of clustering function. "single", "complete" or "average". 
#' @param dist.type Type of distance function. Choices are "spearman", "pearson", or "euclid". 
#'
#' @return An m vector of cluster assignments of the m ordered variables. 
#' @export
#'
#' @examples
#' 
#' data(betas.7)
#' data(annot.7)
#' dat.7.ord <- order.betas.by.chrom.location(betas.7, annot = annot.7)
#' cluster.vec <- Acluster(ordr.vec = dat.7.ord$betas.by.chrom[[1]], thresh.dist = 0.2, 
#' location.vec = dat.7.ord$sites.locations.by.chrom[[1]]$Coordinate_37, max.dist = 1000, type = "average")
#' cluster.vec[1:10]
#' 
Acluster <-
function(ordr.vec, thresh.dist, which.clust = NULL, location.vec = NULL, max.dist = Inf, type = "single", dist.type = "spearman"){
	stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman", "pearson","euclid" )))
	## location.vec gives genomic locations of the sites on ordr.vec.
	## max.dist gives the maximal distance of two adjacent sites able to cluster together
	
	if (!is.null(max.dist) & !all(!is.na(location.vec))) stop("missing location values")
	
	
	### initialize variables
	le <- dim(ordr.vec)[2]
	dist.clust <- rep(0, le - 1)

	if (is.null(which.clust)) { ## create initial clustering indicators vector 
		which.clust <- 1:le
	
	### calculate the first distances vector
	dist.clust <- calc.dist.d.neighbor(ordr.vec, 1, dist.type)

	} else{ ## calculate distances between clusters 
		
		dist.clust <- rep(0, le -1)
		
		for (k in 1:(le - 1)){
			if (which.clust[k] == which.clust[k + 1]) dist.clust[k] <- 0 else{
				clust.1.min <- first(which.clust[which.clust == which.clust[k]])
				clust.1.max <- k
				
				clust.2.min <- k +1
				clust.2.max <- last(which.clust[which.clust == which.clust[k + 1]])
				
				dist.clust[k] <- calc.dist.clusters(ordr.vec[,clust.1.min:clust.1.max], 
                                            ordr.vec[,clust.2.min:clust.2.max], 
                                            type = type, dist.type = dist.type)
				
				}
		}
		}		
	
	### iterate - merging clusters until no more to merge
	no.more <- FALSE
	
	condition1 <- function(x){
		return((x < thresh.dist) & is.null(max.dist) ) 
	}
	
	condition2 <- function(x, loc1, loc2){
		if (is.null(max.dist)) return(FALSE) 
		if (is.null(loc1) | is.null(loc2)) return(FALSE)
		else return((x < thresh.dist) & !is.null(max.dist) & (max.dist >= abs(loc2 - loc1)))
	}
	
	while (!no.more){
		no.more = T
		clust.1.min <- 1
		for (k in 1:(le - 1)){
			if (which.clust[k] != which.clust[k +1 ]){
			if (condition1(dist.clust[k]) | condition2(dist.clust[k], location.vec[k], location.vec[k+1])){ #(dist.clust[k] < thresh.dist )
				no.more = F ## we now merge, and potentially there will be more merging to do
				
				### merge clusters:
				clust.to.merge <- which.clust[k + 1]
					 l <- k+1
					 while (which.clust[l] == clust.to.merge){
						which.clust[l] <- which.clust[k]
						l <- l+1
						if (l > le) break
					}
					clust.1.max <- l - 1
				dist.clust[clust.1.min:(clust.1.max -1)] <- Inf ## so we won't try to merge within cluster	
					 
				
				##  recalculate distances with the following cluster:
				if (clust.1.max < le) {  ## only if there is another cluster
					clust.2.min <- l
					clust.2.max <- l
				
					if (l < le) { 
						l <- l + 1
						  while (which.clust[l] == which.clust[clust.2.min]){
							l <- l+1
							if (l > le) break
						}
						clust.2.max <- l - 1 }
				
					dist.clust[clust.1.max] <- calc.dist.clusters(as.matrix(ordr.vec[,clust.1.min:clust.1.max]), as.matrix(ordr.vec[,clust.2.min:clust.2.max]), type = type, dist.type = dist.type) }
					
				## recalcualte distances with the previous cluster: 
				if (clust.1.min > 1){
					clust.0.max <- clust.1.min - 1
					clust.0.min <- clust.0.max
					
					if (clust.0.min > 1){ l <- clust.0.min - 1
					
					 while (which.clust[l] == which.clust[clust.0.max]){
						which.clust[l] <- which.clust[clust.0.max]
						l <- l-1
						if (l < 1) break
					}
					clust.2.min <- l + 1  }
					
					dist.clust[clust.0.max] <- calc.dist.clusters(as.matrix(ordr.vec[,clust.0.min:clust.0.max]), as.matrix(ordr.vec[,clust.1.min:clust.1.max]), type = type, dist.type = dist.type) }		
				
			}  else clust.1.min <- k + 1 
			}}	
		}
	
	return(which.clust)
}
