#' Title Annotate vector of Illumina probe names
#'
#' Returns a matrix with annotations of the CpG sites represented by a given vector of probe names
#'
#' @param probe.vec A vector of Illumina probe names
#' @param annot A preloaded data.table of annotation
#' @param annotation.file.name If annot is not given, annotation.file.name can provide the file name of the illumine annotation file and the function will load and convert it to a data.table. 
#' If neither annot nor annotation.file.name are provided, annotation will be loaded from Tim Triche's bioconductor R package. 
#' @param required.annotation Names of columns from Illumina annotation file, to be provided as annotation 
#'
#' @return A matrix of requred.annotations for each of the probes in probe.vec
#' @export
#'
#' @examples
#' 
#' data(betas.7)
#' data(annot.7)
#' annot.sites <- annot.probe.vec(rownames(betas.7)[1:10], annot = annot.7)
#' 
#' 
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
