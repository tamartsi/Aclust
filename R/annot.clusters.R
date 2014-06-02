annot.clusters <-
function(clusters.list , annot = NULL, annotation.file.name = NULL, required.annotation = c("IlmnID", "Coordinate_36", "Gene_Name","UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")){
## function that getes a list of clusters, and returns a list of annotations

	require(data.table)
	if (is.null(annot)) {
		if (!is.null(annotation.file.name)){
			cat("Loading annotation from Illumina's menifest", "\n")
			annot <- read.csv(annotation.file.name, skip = 7)
			annot <- data.table(annot)
			setkeyv(annot, c("CHR","Coordinate_36") ) } else{
				cat("Loading annotation from Tim Triche's package on Bioconductor", "\n")
				annot <- create.annot.triche(probe.vec = unlist(clusters.list))
				}
			}

	l <- length(clusters.list)
	clusters.annot.list <- vector(mode = "list", length = l)

	clusters.annot.list <- lapply(clusters.list, function(x) annot.probe.vec(x, annot = annot, annotation.file.name = annotation.file.name, required.annotation = required.annotation))

	return(clusters.annot.list)
}
