print.site.analysis <-
function(res.mat, file.name){
## res.mat is  a matrix with individual sites results, of the sites in the top clusters. 
## The table summarizing the results will be printed to the file file.name	
	
	sink(file.name, append = T)
	cat("\\begin{table}[htbp]", "\n")
	cat("\\begin{center}","\n")
	cat("\\begin{tabular}{", rep("l", ncol(res.mat)), "}", "\n")
	cat("\\hline\\hline","\n")
	
	str.temp <- "\\multicolumn{1}{c}{{\\bf"
	for (i in 1:ncol(res.mat)){
		temp <- gsub("_"," ", colnames(res.mat)[[i]])
		str.temp <- paste(str.temp, temp, "}} ")
		 
		 if (i < ncol(res.mat)) str.temp <- paste(str.temp, "&\\multicolumn{1}{c}{{\\bf") else
		 		str.temp <- paste(str.temp, "\\tabularnewline")
		
	}
	
	cat(str.temp, "\n")
	cat("\\hline", "\n")
	
	
	for (i in 1:nrow(res.mat)){
		temp <- as.character(res.mat[i,1])
		
		for (j in 2:ncol(res.mat)){
			if (is.numeric(res.mat[i,j])) temp <- paste(temp, "&", round(res.mat[i,j],3))
				else temp <- paste(temp, "&", res.mat[i,j])	
		}
		temp <- paste(temp, "\\\\")
		cat(temp, "\n")
		
	}
	

	
	cat("\\hline", "\n")
	
	cat("\\end{tabular}","\n")
	cat("\\caption{Individual site analysis results, for sites in the", length(unique(res.mat[,grep("cluster", colnames(res.mat))])), "top clusters associated with the exposure. The \\textit{P}-value was obtained through the sandwich estimator of the variance}","\n")
	cat("\\end{center}","\n")
	cat("\\end{table}","\n")

	sink()	

}
