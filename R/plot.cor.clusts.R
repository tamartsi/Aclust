plot.cor.clusts <-
function(cor.mat, clusts, erase.non.clust = F, lwd = 2, main = "", xlab = "", labels = NULL){
	## cor.mat is a correlation matrix
	## clusts is a vector indicating for which cluster each variables belong. 
	## assumption: only neighboring variables could be clustered. 
	
	
	square.lines <- function(bottom.left, edge.length){
	## produces the edges of the square for the ploting function 'segments"
	
	from.x <- c(bottom.left, bottom.left, bottom.left + edge.length, bottom.left + edge.length)
	from.y <- c(bottom.left, bottom.left, bottom.left + edge.length, bottom.left + edge.length)
	to.x   <- c(bottom.left, bottom.left + edge.length, bottom.left, bottom.left + edge.length)
	to.y   <- c(bottom.left + edge.length, bottom.left, bottom.left + edge.length, bottom.left)
	
	return(list(x0 = from.x, y0 = from.y, x1 = to.x, y1 = to.y))
}
	
	
	
	if(clusts[1] >1) clusts <- clusts - clusts[1] + 1
	
	if (clusts[2] - clusts[1] > 1) clusts[-1] <- clusts[-1] - (clusts[2] - clusts[1] -1)
	
	clust.sizes <- table(clusts)
	clust.sizes <- clust.sizes[order(clust.sizes, decreasing = T)]
#	n.clust <- min(which(clust.sizes <= 1)) - 1
	n.clust <- length(which(clust.sizes > 1))
	
	if (n.clust == 0 ) {
			 if (is.null(labels)) image(cor.mat, main= main, xlab = xlab) else{
			 	image(cor.mat,  axes = F)
				axis(1, at = seq(0,1,length =  length(labels)), labels = labels, xlab = xlab)
				axis(2, at = seq(0,1,length =  length(labels)), labels = labels, xlab = xlab)
			 	
			 }
		return()}
	
	## parameter for the image
	sq <- 1/(dim(cor.mat)[1] - 1)
	half.sq <- sq/2
	qrt.sq <- sq/4
	
	cor.clus.mat <- cor.mat
	segment.list <- vector("list", length = n.clust)
	ind.clus <- 1
	ind.line <- 1
	 while (ind.clus <= dim(cor.mat)[1]){
		clust <- which(clusts == ind.clus)
		if (length(clust) == 1){
			if (erase.non.clust) cor.clus.mat[ind.clus, - ind.clus ] <- cor.clus.mat[- ind.clus, ind.clus] <- NA
			ind.clus <- ind.clus + 1
			} else {
				lc <- length(clust)
				bottom.left <- (ind.clus-1)*sq - sq/2
				edge.length <- lc*sq
				segment.list[[ind.line]] <- square.lines(bottom.left, edge.length)
				ind.line <- ind.line + 1
				
				if (erase.non.clust) cor.clus.mat[clust, -clust] <- cor.clus.mat[-clust, clust] <- NA
			
				ind.clus <- ind.clus + lc
			}
			
		}
		
		if (!is.null(main)){
			if (is.null(labels)) image(cor.clus.mat, main = main, xlab = xlab) else{
				image(cor.clus.mat, main = main, axes = F, xlab = xlab)
				axis(1, at = seq(0,1,length = length(labels)), labels = labels)
				axis(2, at = seq(0,1,length =  length(labels)), labels = labels)
				}
		}  else{
			 if (is.null(labels)) image(cor.clus.mat, xlab = xlab) else{
			 	image(cor.clus.mat,  axes = F, xlab = xlab)
				axis(1, at = seq(0,1,length =  length(labels)), labels = labels)
				axis(2, at = seq(0,1,length =  length(labels)), labels = labels)
			 	
			 }
			}
			
		for (i in 1:length(segment.list)){
			segments(segment.list[[i]]$x0, segment.list[[i]]$y0, segment.list[[i]]$x1, segment.list[[i]]$y1, lwd = lwd)
		}
		
	}
