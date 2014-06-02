print.clusters.gee.summary <-
function(clusters.res, annotation.clusters, file.name){
	## clusters.res is a matrix that summrizes the exposure effect on nrow(clusters.res) clusters. 
	## annotation matrix has the annotation of the clusters
	sink(file.name, append = T)
	cat("\\begin{sidewaystable}[htbp]", "\n")
		if (ncol(annotation.clusters) > 6) cat("\\small", "\n")
	cat("\\begin{center}","\n")
	cat("\\begin{tabular}{", rep("l", ncol(annotation.clusters)), "}", "\n")
	cat("\\hline\\hline","\n")
	str.temp <- "\\multicolumn{1}{c}{{\\bf"
	for (i in 1:ncol(annotation.clusters)){
		temp <- gsub("_"," ", colnames(annotation.clusters)[[i]])
		str.temp <- paste(str.temp, temp, "}} ")
		 
		 if (i < ncol(annotation.clusters)) str.temp <- paste(str.temp, "&\\multicolumn{1}{c}{{\\bf") else
		 		str.temp <- paste(str.temp, "\\tabularnewline")
		
	}
	
	cat(str.temp, "\n")
	cat("\\hline", "\n")
	
	
	## now print each cluster
	for (i in 1:nrow(clusters.res)){
		temp.str <- " \\multicolumn{2}{c}{Exposure effect: "
		temp.str <- paste(temp.str, round(clusters.res$exposure_effect_size[i],2), "}& \\multicolumn{2}{c}{Exposure {\\it P}-value: ")
			#and.vec <- rep("&" , ncol(annotation.clusters) - 4)
			and.vec <- "&"
			if (ncol(annotation.clusters) - 5 > 0){
				for (k in 1:(ncol(annotation.clusters) - 5)) and.vec <- paste(and.vec, "&")			
			}
			
		temp.str <- paste(temp.str, round(clusters.res$fdr_pval[i],3), "}", and.vec, "\\\\ \\hline"  )
		cat(temp.str, "\n")
		
		ind.sites <- which(annotation.clusters$cluster == i) 
		## print every site in the cluster:
		for (j in ind.sites){
			temp <- as.character(annotation.clusters[j,1])
			for (k in 2:ncol(annotation.clusters)){
				if (k < ncol(annotation.clusters)) temp <- paste( temp, "&", gsub("_", " ", annotation.clusters[j,k])) else
					temp <- paste(temp, "&", gsub("_", " ", annotation.clusters[j,k]), "\\\\ ")
					sub("_", " ", temp)
			}
			cat(temp, "\n")	
		}
		
		cat("\\hline", "\n")
	
	}
	
	cat("\\hline", "\n")
	
	
	cat("\\end{tabular}","\n")
	cat("\\caption{The sites in the", nrow(clusters.res), "top clusters associated with the exposure selected in the analysis.}","\n")
	cat("\\end{center}","\n")
	cat("\\end{sidewaystable}","\n")
	if (ncol(annotation.clusters) > 6) cat("\\normalsize", "\n")
	sink()	
	
}
