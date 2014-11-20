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
