#' Title Load annotation table 
#' 
#' The package is adjusted to Illumina annotation .csv file. 
#' 
#' @param annotation.file.name The name of the .csv Illumina annotation file. 
#'
#' @return A data.table object with the annotation.
#' @export
#'
#' @examples
#' 
#' annotation.file.name <- "illumina_450_manifest_v.1.2.csv"
#' annot <- load.annotation(annotation.file.name)
#' annot[1:5]
#' 
load.annotation <-
function(annotation.file.name){
	annot <- read.csv(annotation.file.name, skip = 7, as.is = T)
	annot <- data.table(annot)
  setnames(annot, c("MAPINFO"), c("Coordinate_37"))
  annot$Coordinate_37 <- as.numeric(as.character(annot$Coordinate_37))
	setkeyv(annot, c("CHR","Coordinate_37"))
	return(annot)
}
