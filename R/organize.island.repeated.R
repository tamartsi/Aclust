organize.island.repeated <-
function(X, covar.mat){
	## function recives X matrix of betas values (columns are people) and matrix covar.mat of covariates 
	## (columns are covariates). returns reshaped matrix. 
	
	data <- cbind(t(X), covar.mat)
	data.long <- reshape(data, direction = "long", varying = list(names(data)[1:dim(X)[1]]), v.names = "Beta", idvar = colnames(covar.mat), timevar="probeID", times = rownames(X))
	rownames(data.long) <- unlist(lapply(rownames(rownames(X)), function(x) rep(x,ncol(X))))    
	
	return(data.long)
	
}
