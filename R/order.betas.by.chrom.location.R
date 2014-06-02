order.betas.by.chrom.location <-
function(betas, annot = NULL, annotation.file.name = NULL, return.chroms = NULL)  {
### function that gets a matrix of beta values (rows are CpG sites, column names are samples' names) 
### and returns a list of matrices ordered for each chromosome, or for selected chromosomes
### also returns a vector of locations for each probe
	require(data.table)
	
	if (is.null(annot)) {
		if (!is.null(annotation.file.name)){
			cat("Loading annotation from Illumina's menifest", "\n")
			annot <- read.csv(annotation.file.name, skip = 7)
			annot <- data.table(annot)
			setkeyv(annot, c("CHR","Coordinate_36") ) } else{
				cat("Loading annotation from Tim Triche's package on Bioconductor", "\n")
				annot <- create.annot.triche(probe.vec = rownames(betas), only.locations = T)
				}
			}
	
	annot.betas <- annot[IlmnID %in% rownames(betas)]
	annot.betas$IlmnID <- as.character(annot.betas$IlmnID)
	annot.betas$Coordinate_36 <- annot.betas$Coordinate_36
	
	chroms <- unique(annot.betas$CHR)
	if (!is.null(return.chroms)) chroms <- intersect(chroms, return.chroms)
	chroms <- chroms[order(as.numeric(as.character(chroms)))]
	betas.by.chrom <- vector(mode = "list", length = length(chroms))
	sites.by.chrom <- vector(mode = "list", length = length(chroms))
	names(betas.by.chrom) <- names(sites.by.chrom) <- chroms
	
	
	for (i in 1:length(chroms)){
		cpg.chrom <- as.character(annot.betas[CHR == chroms[i]]$IlmnID)
		betas.by.chrom[[i]] <- as.matrix(betas[cpg.chrom,])
		if (ncol(betas.by.chrom[[i]]) == 1){
			betas.by.chrom[[i]] <- t(betas.by.chrom[[i]])
			rownames(betas.by.chrom[[i]]) <- cpg.chrom
		}

		sites.by.chrom[[i]] <- annot.betas[CHR == chroms[i], c("IlmnID", "Coordinate_36"), with = F]
		
	}
	
	return(list(betas.by.chrom = betas.by.chrom, sites.locations.by.chrom = sites.by.chrom))
}
