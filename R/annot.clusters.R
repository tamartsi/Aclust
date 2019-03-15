#' Title Annotate vectors of Illumina probe names, each vector represent different cluster
#' 
#' Returns a list of matrices with annotations of the CpG sites represented by vector of probe names. Each vector represent a cluster
#'
#' @param clusters.list A list in which each element is a vector of illumine probe names. 
#' @param annot A preloaded data.table of annotation
#' @param annotation.file.name If annot is not given, annotation.file.name can provide the file name of the illumine annotation file and the function will load and convert it to a data.table. 
#' If neither annot nor annotation.file.name are provided, annotation will be loaded from Tim Triche's bioconductor R package. 
#' @param required.annotation Names of columns from Illumina annotation file, to be provided as annotation 
#'
#' @return A list corresponding to clusters.list. Each element in the list is matrix of required.annotations  for each of the probes in the vector in the corresponding element in clusters list
#' @export
#'
#' @examples
#' 
#' data(betas.7) ## upload methylation data
#' exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
#' covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
#' data(annot.7)
#' clusters.list <- assign.to.clusters(betas.7, annot.7)
#' annotated.clusters <- annot.clusters(clusters.list[1:10], annot = annot.7, required.annotation = c("IlmnID", "Coordinate_37", "UCSC_RefGene_Name", "UCSC_RefGene_Group"))
#' 

annot.clusters <-
function(clusters.list, annot = NULL, 
         annotation.file.name = NULL, 
         required.annotation = c("IlmnID", "Coordinate_37", 
                                 "UCSC_RefGene_Name","UCSC_RefGene_Group", 
                                 "UCSC_CpG_Islands_Name", 
                                 "Relation_to_UCSC_CpG_Island")){
## function that gets a list of clusters, and returns a list of annotations

	if (is.null(annot)) {
		if (!is.null(annotation.file.name)){
			cat("Loading annotation from Illumina's manifest", "\n")
			annot <- read.csv(annotation.file.name, skip = 7)
			annot <- data.table(annot)
      setnames(annot, c("MAPINFO"), c("Coordinate_37"))
			setkeyv(annot, c("CHR","Coordinate_37") ) } else{
				cat("Loading annotation from Tim Triche's package on Bioconductor", "\n")
				annot <- create.annot.triche(probe.vec = unlist(clusters.list))
				}
			}

	l <- length(clusters.list)
	clusters.annot.list <- vector(mode = "list", length = l)

	clusters.annot.list <- lapply(clusters.list, function(x) annot.probe.vec(x, annot = annot, annotation.file.name = annotation.file.name, required.annotation = required.annotation))

	return(clusters.annot.list)
}
