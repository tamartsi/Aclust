#' Title GEE analysis of exposure effect on clusters of CpG sites
#'
#' The function takes a list of clusters (given as vector of probe names), and perform a GEE analysis of the beta values corresponding to these clusters. 
#' i.e. for each cluster of probes a GEE model is fit, were the sites in the clusters are treated as multiple outcomes. 
#' The analysis assumes a common exposure effect on all outcomes, but different intercept for each sites. 
#'
#' @param betas An (m by n) matrix of methylation values of n participants measured on m methylation sites. 
#' @param clusters.list A list where each item is a vector of methylation sites, composing a single cluster.
#' @param exposure The exposure of interest, a vector ordered according to the columns of the matrix betas. 
#' @param covariates An (n by p) matrix of covariates. The order of rows should match the order of columns of the matrix betas. 
#' @param id A vector of ids of $n$ participants. The order should match the order of columns of the matrix betas. 
#' @param working.cor Assume working correlation structure. Could be "unstructured" (not recommended) "independence", "exchangeable", or "ar1".
#' @param minimum.cluster.size ???
#' @param result.file.name Name of file to print analysis results to
#'
#' @return data frame where each row provides a vector with probes in the cluster, the estimated effect size of exposure on the cluster, estimated (sandwich) standard error, p-value, number of sites in the cluster, and the chromosome where the cluster is.
#' @export
#'
#' @examples
#' 
#' 
#' data(betas.7) ## upload methylation data
#' exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
#' covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
#' rownames(covariates) <- colnames(betas.7)
#' data(annot.7)
#' clusters.list <- assign.to.clusters(betas.7, annot.7)
#' GEE.results.clusters <- GEE.clusters(betas.7, clusters.list, exposure, covariates, id = colnames(betas.7), working.cor = "ex")
#' top.clusters.summary <- summarize.top.clusters(betas.7, covariates, exposure, id = colnames(betas.7), GEE.results.clusters, "results.tex", annot= annot.7)
#' 
GEE.clusters <-
function(betas, clusters.list, exposure, covariates, id, 
         working.cor = "ex", minimum.cluster.size = 2, 
         result.file.name = NULL){
## the ids are rownames of the betas, and are ordered by the same ordering of the covariates matrix.
## function that gets a matrix of beta values, a list of clusters, exposure variable, and covariates
## and performs GEE. Returns a matrix with the results of the GEE model (only the exposure effect and p-value)
## for each cluster. If the cluster has a single probe?
## The results are ordered according to p-value
## Note: methylation is treated here as outcome. 
## if file is given, print results to file. 

	n.sites <- unlist(lapply(clusters.list, function(x) return(length(x))))
	inds.rm <- which(n.sites < minimum.cluster.size)
	
	if (length(inds.rm) > 0) clusters.list <- clusters.list[-inds.rm]
	
	
	n.mod <- length(clusters.list)
	
	effect <- rep(0, n.mod)
	se <- rep(0, n.mod)
	pvals <- rep(0, n.mod)
	sites <- rep("", n.mod)
	n.sites <- rep(0, n.mod)

	if (is.null(dim(covariates))) covariates <- as.matrix(covariates)

	if (is.null(colnames(covariates)))  colnames(covariates) <- rep("", ncol(covariates))
	for (i in 1:ncol(covariates)){
		if (colnames(covariates)[i] == "") colnames(covariates)[i] <- paste("covariate", i, sep = '_')
		
	}
	
	model.expression <- "model <- geeglm(Beta ~ exposure +"
	for (j in 1:ncol(covariates)){
		 model.expression <- paste(model.expression, colnames(covariates)[j], "+") 	}

	model.expression <- paste(model.expression, "as.factor(probeID), id = as.numeric(id), data = temp.long, corstr = working.cor)")
	
	model.expr.1.site <- paste("model <- geeglm(clus.betas[ind.comp] ~ exposure[ind.comp] + ")	
	for (j in 1:ncol(covariates)){
		 if (j < ncol(covariates)) model.expr.1.site <- paste(model.expr.1.site, colnames(covariates)[j], "+")
		 	else  	model.expr.1.site <- paste(model.expr.1.site, colnames(covariates)[j])}
		 	
	model.expr.1.site <- paste(model.expr.1.site, ", data = covariates[ind.comp,], id = as.numeric(id)[ind.comp])")
	covariates <- as.data.frame(covariates)

	for (i in 1:n.mod){
		clus.probes <- clusters.list[[i]]
		clus.betas <- betas[clus.probes,]
			
			
		if (length(clus.probes) == 1){
			
			ind.comp <- which(complete.cases(cbind(clus.betas, exposure, covariates, id)))
			eval(parse(text = model.expr.1.site))
			
			effect[i] <- summary(model)[[6]][2,1]
			se[i] <- summary(model)[[6]][2,2]
			pvals[i] <- summary(model)[[6]][2,4]
			n.sites[i] <- 1
			sites[i] <- clus.probes
			
			
		}	else{
			
		ind.comp <- which(complete.cases(cbind(t(clus.betas), exposure, covariates, id)))
			
		temp.betas <- clus.betas[,ind.comp]
		
		temp.covars <- cbind(exposure[ind.comp], covariates[ind.comp,], id[ind.comp])
		colnames(temp.covars) <- c("exposure", colnames(covariates), "id")
			
		temp.long <- organize.island.repeated(temp.betas, temp.covars)
		temp.long <- temp.long[complete.cases(temp.long),]
		temp.long <- temp.long[order(temp.long$id),]

		eval(parse(text = model.expression))
				
		
		effect[i] <- summary(model)[[6]][2,1]
		se[i] <- summary(model)[[6]][2,2]
		pvals[i] <- summary(model)[[6]][2,4]
		n.sites[i] <- length(clus.probes)
		for (c in 1:length(clus.probes)) sites[i] <- paste(sites[i], clus.probes[c], sep = ";")
		substr(sites[i], 1, 1) <- ""
		
		}
		
	}

	result <- data.frame(exposure_effect_size = effect, 
                       exposure_effect_se = se, 
                       exposure_pvalue= pvals, 
                       n_sites_in_cluster = n.sites, 
                       cluster_sites =sites)	
	result <- result[order(result$exposure_pvalue),]
	
	
	if(!is.null(result.file.name)){
    write.table(result, file = result.file.name, 
                col.names = T , row.names = F,  append = T)
	}	
	return(result)

}
