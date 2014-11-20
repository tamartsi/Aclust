annot.probe.vec <-
function(probe.vec, annot = NULL, 
         annotation.file.name = NULL, 
         required.annotation = c("IlmnID", "Coordinate_37", 
                                 "UCSC_RefGene_Name","UCSC_RefGene_Group", 
                                 "UCSC_CpG_Islands_Name", 
                                 "Relation_to_UCSC_CpG_Island")){
## function that gets a vector of Illumina ids and returns
## a data.frame of annotation according to the required annotation 
## one can specify column names from the illumina annotation file). 

  if (is.null(annot)) {
		if (!is.null(annotation.file.name)){
			cat("Loading annotation from Illumina's manifest", "\n")
			annot <- read.csv(annotation.file.name, skip = 7)
			annot <- data.table(annot)
      setnames(annot, c("MAPINFO"), c("Coordinate_37"))
			setkeyv(annot, c("CHR","Coordinate_37") ) } else{
				cat("Loading annotation from Tim Triche's package on Bioconductor", "\n")
				annot <- create.annot.triche(probe.vec = probe.vec)
				}
			}
	
	t1 <- grep("UCSC_RefGene_Name", required.annotation)		
	t2 <- grep("UCSC_RefGene_Name", names(annot))
	if (length(t1) > 0) required.annotation[t1] <- names(annot)[t2]
	
	annot.probes <- annot[match(probe.vec, annot$IlmnID), required.annotation, with = F]

	return(as.data.frame(annot.probes))
}
