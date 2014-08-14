annot.probe.vec <-
function(probe.vec, annot = NULL, annotation.file.name = NULL, required.annotation = c("IlmnID", "Coordinate_36", "Gene_Name","UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")){
## function that gets a vector of Illumina ids and returs a matrix of annotations according to the required annotation 
## one can specify column names from the illumina annotation file). 

  if (is.null(annot)) {
		if (!is.null(annotation.file.name)){
			cat("Loading annotation from Illumina's menifest", "\n")
			annot <- read.csv(annotation.file.name, skip = 7)
			annot <- data.table(annot)
			setkeyv(annot, c("CHR","Coordinate_36") ) } else{
				cat("Loading annotation from Tim Triche's package on Bioconductor", "\n")
				annot <- create.annot.triche(probe.vec = probe.vec)
				}
			}
	
	t1 <- grep("Gene_Name", required.annotation)		
	t2 <- grep("Gene_Name", names(annot))
	if (length(t1) > 0) required.annotation[t1] <- names(annot)[t2]
	
	annot.probes <- annot[match(probe.vec, annot$IlmnID), required.annotation, with = F]

	return(as.data.frame(annot.probes))

}
