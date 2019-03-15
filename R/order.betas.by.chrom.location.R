#' Title Organizes betas by location
#' 
#' Organizes a matrix of methylation values to be ordered by chromosomal location
#' 
#' @param betas An (m by n) matrix of methylation values measured on $n$ participants in $m$ sites. 
#' @param annot A preloaded data.table of annotation
#' @param annotation.file.name If annot is not given, annotation.file.name can provide the file name of the illumine annotation file and the function will load and convert it to a data.table. 
#' If neither annot nor annotation.file.name are provided, annotation will be loaded from Tim Triche's bioconductor R package.
#' @param return.chroms Optional list of chromosomes, if one is interested in specific chromosomes. 
#'
#' @return
#' \item{betas.by.chrom}{A list ordered by chromosome number. Each item in this list contains a matrix of methylation values for the subsets of sites from betas from the corresponding chromosome. The rows are organized according to chromosomal location. }
#' \item{sites.locations.by.chrom}{A list ordered by chromosome number. Each item in this list contains a matrix specifying chromosomal locations of each of the sites in the subsets of sites from betas from the corresponding chromosome. The rows are organized according to chromosomal location. }
#'  
#' @export
#'
#' @examples
#' 
#' data(betas.7) ## upload methylation data
#' data(annot.7)
#' dat.7.ord <- order.betas.by.chrom.location(betas.7, annot = annot.7)
#' dat.7.ord$betas.by.chrom$"7"[1:5,1:5]
#' dat.7.ord$sites.locations.by.chrom$"7"[1:5]
#' 
order.betas.by.chrom.location <-
function(betas, annot = NULL, annotation.file.name = NULL, return.chroms = NULL)  {
### function that gets a matrix of beta values (rows are CpG sites, column names are samples' names) 
### and returns a list of matrices ordered for each chromosome, or for selected chromosomes
### also returns a vector of locations for each probe
	
  if (is.null(annot)) {
		if (!is.null(annotation.file.name)){
			cat("Loading annotation from Illumina's manifest", "\n")
			annot <- read.csv(annotation.file.name, skip = 7)
			annot <- data.table(annot)
			setnames(annot, c("MAPINFO"), c("Coordinate_37"))
      setkeyv(annot, c("CHR","Coordinate_37") ) } else{
				cat("Loading annotation from Tim Triche's package on Bioconductor", "\n")
				annot <- create.annot.triche(probe.vec = rownames(betas), only.locations = T)
				}
			}
	
	annot.betas <- annot[IlmnID %in% rownames(betas)]
	annot.betas$IlmnID <- as.character(annot.betas$IlmnID)
	annot.betas$Coordinate_37 <- annot.betas$Coordinate_37
	
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

		sites.by.chrom[[i]] <- annot.betas[CHR == chroms[i], 
                                       c("IlmnID", "Coordinate_37"), 
                                       with = F]
		
	}
	
	return(list(betas.by.chrom = betas.by.chrom, 
              sites.locations.by.chrom = sites.by.chrom))
}


#Returns A list with two components: 1) Matrices of methylation values for the subsets of sites from betas from the corresponding chromosome. The rows are organized according to chromosomal location.
#2) Matrices specifying chromosomal locations of each of the sites in the subsets of sites from betas from the corresponding chromosome. The rows are organized according to chromosomal location.
