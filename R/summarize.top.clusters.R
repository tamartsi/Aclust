#' Title Summarize the top clusters
#' 
#' Returns information about the top cluster in the GEE analysis
#' 
#' @param betas An (m by n) matrix of methylation values of $m$ methylation sites measured on $n$ individuals
#' @param covariates Either a (n by p) vectors of adjusting covariates, ordered so that its rows corresponds to the columns of the `betas' matrix, or NULL if there are no covariates.  
#' @param exposure A vector of size $n$ of exposure values for each individual. Same order as the covariates. 
#' @param id A vector of size $n$ of the IDs of the individual. The IDs should match the column names of the matrix `betas'.  
#' @param clusters.GEE.results A matrix of results from a  GEE analysis of clusters.
#' @param clusters.GEE.results.file A file with the results from a GEE analysis of clusters (if clusters.GEE.results is not given)
#' @param minimum.sites A minimum number of sites in a cluster to be considered in the analysis. 
#' @param top.number The required number of top clusters to be reported. 
#' @param cutoff.fdr.pval A significance p-value threshold (after FDR correction). Only clusters with p-value smaller than cutoff.fdr.pval will be reported. 
#' @param cutoff.effect.size An effect size threshold. Only clusters with estimated effect size larger than cutoff.effect.size will be reported. 
#' @param annot annotation data table. The package uses the Illumina annotation file ``illumina_450_manifest_v.1.2". 
#' @param annotation.file.name A name of annotation file to read. By default it is not required, since one can use `annot', which is in the package. 
#' @param required.annotation What is the annotation to be reported on each of the sites in the reported clusters? 
#' @param file.to.print.report File name to print tables of annotation, effect sizes, and individual site analyses of sites in the top clusters. 
#' @param print.progress Print status messages while the function progresses
#'
#' @return
#' \item{top.clusters}{A data frame summarizing the GEE results (exposure effect estimates, etc) of the top clusters }
#' \item{annotation.top.clusters}{Annotations of the sites from the top clusters}
#' \item{individual.sites.analysis}{Results of individual sites analysis, of the sites in the top clusters.  }

#' @export
#'
#' @examples
summarize.top.clusters <-
function(betas, covariates, exposure, id, clusters.GEE.results = NULL, clusters.GEE.results.file = NULL, minimum.sites = 2, top.number = 10, cutoff.fdr.pval = 0.05, cutoff.effect.size = NULL, annot = NULL, annotation.file.name = NULL, required.annotation = c("IlmnID", "Coordinate_37", "UCSC_RefGene_Name","UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"), file.to.print.report = NULL, print.progress = F){
## Returns a summary of the top clusters:
## Chooses the top clusters according to the top.number provvided, with FDR p-value at most 0.05, and with effect size at least cutoff.effect.size (if given)
## For each of the top cluster, returns exposure effect size and pvalue before and after FDR, and illumina annotation
## Also returns in dividual GEE analysis of sites in the cluster (i.e. with sandwich standard error)
## Finally, the returnd lists are ordered by exposure effect size on the cluster
## If a file to print report is given, then a latex file is generated with two tables. One table summarized the annotation and effects of
## the top clusters, annother table summarizes the individual sites analysis of the sites from the top clusters. 
## minimum.sites is the minimal size of clusters to consider. The FDR correction is applied only these clusters. 
	require(data.table)
	 
	if (!is.null(clusters.GEE.results)) anal.results <- clusters.GEE.results 
		else anal.results <- read.table(clusters.GEE.results.file, header = T, as.is = T)
	
			
	anal.results$cluster_sites <- as.character(anal.results$cluster_sites)	
	anal.results <- data.table(anal.results)
	anal.results <- anal.results[n_sites_in_cluster >= minimum.sites]



	anal.results[,fdr_pval := p.adjust(anal.results$exposure_pvalue, "fdr")]

#	anal.results[,fdr_pval := 0.04]  ## to test the function
	
	anal.results[,abs.exp.effect := abs(exposure_effect_size)]
	setkeyv(anal.results, cols = "abs.exp.effect")
	
	anal.results <- anal.results[fdr_pval <= cutoff.fdr.pval]
	if (!is.null(cutoff.effect.size)) anal.results <- anal.results[exposure_effect_size >= cutoff.effect.size]
	
	top.number <- min(top.number, nrow(anal.results))	
	
	if (print.progress) message(paste("up to ", top.number, "of clusters will be analyzed and reported"))
	
	if (top.number == 0){
		cat("0 clusters were found with corrected p-value smaller than ",  cutoff.fdr.pval, "no output produced" ,"/n")
	}  else {
	
	
	
	last.one <- nrow(anal.results)	
	top.clusters <- anal.results[last.one:(last.one - top.number + 1)]
	
	### Top.clusters will be outputed 
	
	
	sites.vec <- rep("", length = sum(top.clusters$n_sites_in_cluster))
	cluster.vec <- rep(0, sum(top.clusters$n_sites_in_cluster))
	
	ind <- 1
	for (i in 1:top.number){
		
		sites.vec[ind:(ind + top.clusters$n_sites_in_cluster[i] -1)] <- strsplit(top.clusters$cluster_sites[i], split = ";")[[1]]
		cluster.vec[ind:(ind + top.clusters$n_sites_in_cluster[i] -1)] <- i

		ind <- ind + top.clusters$n_sites_in_cluster[i]
		
		if (print.progress) message(paste("The ", i, "-th cluster sites were parsed"))
		
	}
	
	annotation.top.clusters <- annot.probe.vec(sites.vec, annot = annot, annotation.file.name = annotation.file.name, required.annotation = required.annotation)
	
	if (print.progress){ 
		message(paste("This is the vector with cluster assignments of the top clusters, of length", length(cluster.vec)))
		print(cluster.vec)
		
		message(paste("This is the annotation of top clusters, describing", nrow(annotation.top.clusters), "sites"))
		print(annotation.top.clusters)}
	
	annotation.top.clusters <- cbind(cluster.vec, annotation.top.clusters)
	colnames(annotation.top.clusters)[1] <- "cluster"
	
	#### annotation.top.clusters is one of the required outputs. 
	
	## finally, inidividual site analysis of the sites in the top clusters:
	ind.res.mat <- matrix(0, nrow = length(sites.vec), ncol = 2)
#	colnames(ind.res.mat) <- c("exp_effect", "exp_pval")
	rownames(ind.res.mat) <- sites.vec
	
	
	model.expr.1.site <- paste("model <- geeglm(temp.meth[ind.comp] ~ exposure[ind.comp]")	
	if (!is.null(covariates)){
	for (j in 1:ncol(covariates)){
		 if (j == ncol(covariates)) model.expr.1.site <- paste(model.expr.1.site, "+", colnames(covariates)[j])
		 	else  	model.expr.1.site <- paste(model.expr.1.site,"+", colnames(covariates)[j])} }
	model.expr.1.site <- paste(model.expr.1.site, ", data = as.data.frame(covariates[ind.comp,]), id = as.numeric(id)[ind.comp])")
	
	if (print.progress) message("starting individual site analysis...")
	
	for (i in 1:length(sites.vec)){

		temp.meth <- betas[sites.vec[i],]
		ind.comp <- which(complete.cases(cbind(temp.meth, exposure, covariates, id)))
		
		eval(parse(text = model.expr.1.site))
		
		ind.res.mat[i,1] <-summary(model)[[6]][2,1]
		ind.res.mat[i,2] <- summary(model)[[6]][2,4]
		
		if (print.progress) print(paste("Analyzed the ", i, "-th site out of ", length(sites.vec), "sites"))

	}
	ind.res.mat <- data.frame(ind.res.mat)
	ind.res.mat <- cbind(sites.vec, cluster.vec, ind.res.mat)
	colnames(ind.res.mat) <- c("site", "cluster", "exposure effect", "exposure p-values")
	### This is the third required output.
	
	
	if (!is.null(file.to.print.report)){
		if (print.progress) message("printing report to file...")
		
		sink(file.to.print.report)
		cat("\\documentclass{article}", "\n")
		cat("\\usepackage{multirow}", "\n")
		cat("\\usepackage{rotating}", "\n")
		cat("\\title{ Cluster analysis report}", "\n")
		cat("\\date{", strsplit(date()," ")[[1]][c(1,2,4,6)], "}", "\n")
		cat("\\begin{document}" , "\n")
		cat("\\maketitle", "\n")
		sink()

		print.clusters.gee.summary(top.clusters, annotation.top.clusters, file.to.print.report)
		print.site.analysis(ind.res.mat, file.to.print.report)

		sink(file.to.print.report, append = T)
		cat("\\end{document}", "\n")
		sink()
		if (print.progress) message("finished preparing latex file with summarizing tables")
	}
	
	top.clusters <- top.clusters[,abs.exp.effect := NULL]
	return(list(top.clusters = top.clusters, annotation.top.clusters = annotation.top.clusters, individual.sites.analysis = ind.res.mat))
	}
}
