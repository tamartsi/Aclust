#' Title Cluster a set of probes
#' 
#' Returns a list with clusters (each item is a cluster - a set of probes)
#' 
#' @param betas An (m by n) matrix of methylation values measured on $n$ participants in $m$ sites. 
#' @param annot A preloaded data.table of annotation
#' @param annotation.file.name If annot is not given, annotation.file.name can provide the file name of the illumine annotation file and the function will load and convert it to a data.table. 
#' If neither annot nor annotation.file.name are provided, annotation will be loaded from Tim Triche's bioconductor R package. 
#' @param dist.type Type of similarity distance function "spearman", "pearson" (correlation measures) or "euclid". 
#' @param method Type of clustering function. Choices are "single", "complete", or "average". 
#' @param dist.thresh A similarity distance threshold. Two neighboring clusters are merged to a single cluster if the similarity distance between them is above dist.thresh. 
#' @param bp.merge A distance in chromosomal location. Any set of methylation sites within an interval smaller or equal to bp.dist will be potentially merged, depending on the similarity between sites at the ends of the interval.
#' @param bp.thresh.clust Optional maximum length between neighboring variables permitting to cluster them together. 
#' @param return.chroms ???
#'
#' @return A list in which each item is a vector of Illumina probe ids of methylation sites clustered together.
#' @export
#'
#' @examples
#' 
#' data(betas.7) ## upload methylation data
#' exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
#' covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
#' 
#' @details The package is adjusted to the Illumina annotation (i.e. through part of its column names) but can easily be adjusted to an attestation data table with different column names. 
#' 
assign.to.clusters <-
function(betas, annot = NULL, annotation.file.name = NULL, dist.type = "spearman", method = "single", dist.thresh = 0.5, bp.merge = 0, bp.thresh.clust = Inf, return.chroms = NULL){
### function that gets a matrix of beta values and returns a list of cluster assignments. it uses the
### function order.betas.by.chrom.location to order the data corresponding to each chromosome and 
### then applies the clustering algorithm	
## return the list ordered by the number of sites in the cluster
	
	chrom.list <- order.betas.by.chrom.location(betas, annot = annot, annotation.file.name = annotation.file.name, return.chroms = return.chroms)
	clusters.by.chrom <- vector(mode = "list", length = length(chrom.list[[1]]))
	
	
	for (i in 1:length(chrom.list[[1]])){
		betas.temp <- chrom.list[[1]][[i]]
		locations.temp <- chrom.list[[2]][[i]]
		
		betas.temp <- betas.temp[which(!is.na(locations.temp$Coordinate_37)),]
		locations.temp <- locations.temp[!is.na(Coordinate_37)]
		
		if (!is.null(bp.merge)){
      which.clust <- Dbp.merge(t(betas.temp),
                               thresh.dist = dist.thresh, 
                               bp.thresh.dist = bp.merge, 
                               as.numeric(locations.temp$Coordinate_37),  
                               dist.type = dist.type)
		}	else which.clust <- 1:nrow(locations.temp)
		if (!is.null(bp.thresh.clust)){	
      clust.vec <- Acluster(t(betas.temp), 
                            thresh.dist = dist.thresh, 
                            which.clust = which.clust, 
                            location.vec = as.numeric(locations.temp$Coordinate_37), 
                            max.dist =  bp.thresh.clust, type = method, 
                            dist.type = dist.type)  
    } else 
					clust.vec <- Acluster(t(betas.temp), thresh.dist = NULL, 
                                which.clust = which.clust, 
                                location.vec = NULL, 
                                max.dist =  bp.thresh.clust, 
                                type = method, dist.type = dist.type)

		clusters.by.chrom[[i]] <- lapply(clust.vec, function(x) {
      return(locations.temp$IlmnID[which(clust.vec == x)])
    })
		clusters.by.chrom[[i]] <- clusters.by.chrom[[i]][which(!duplicated(clusters.by.chrom[[i]]))]
		
		if (i == 1) clusters.all <- clusters.by.chrom[[1]] 
			else clusters.all <- c(clusters.all, clusters.by.chrom[[i]])

	}
	
	return(clusters.all)
	
}
