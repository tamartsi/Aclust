load.annotation <-
function(annotation.file.name){
	annot <- read.csv(annotation.file.name, skip = 7, as.is = T)
	annot <- data.table(annot)
  setnames(annot, c("MAPINFO"), c("Coordinate_37"))
  annot$Coordinate_37 <- as.numeric(as.character(annot$Coordinate_37))
	setkeyv(annot, c("CHR","Coordinate_37"))
	return(annot)
}
