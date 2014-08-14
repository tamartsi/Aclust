load.annotation <-
function(annotation.file.name){
	annot <- read.csv(annotation.file.name, skip = 7, as.is = T)
	annot <- data.table(annot)
	annot$Coordinate_36 <- as.numeric(as.character(annot$Coordinate_36))
	setkeyv(annot, c("CHR","Coordinate_36"))
	return(annot)
}
