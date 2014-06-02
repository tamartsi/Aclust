\name{GEE.clusters}
\alias{GEE.clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs GEE analysis of exposure effect on clusters of CpG sites
}
\description{
The function takes a list of clusters (given as vector of probe names), and perform a GEE analysis of the beta values corresponding to these
clusters. i.e. for each cluster of probes a GEE model is fit, were the sites in the clusters are treated as multiple outcomes. The analysis assumes
a common exposure effect on all outcomes, but different intercept for each sites. 
}
\usage{
GEE.clusters(betas, clusters.list, exposure, covariates, id, working.cor = "ex", result.file.name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
An (m by n) matrix of methylation values of n participants measured on m methylation sites. 
}
  \item{clusters.list}{
A list where each item is a vector of methylation sites, composing a single cluster.
}
  \item{exposure}{
The exposure of interest, a vector ordered according to the columns of the matrix betas. 
}
  \item{covariates}{
An (n by p) matrix of covariates. The order of rows should match the order of columns of the matrix betas. 
}
  \item{id}{
A vector of ids of $n$ participants. The order should match the order of columns of the matrix betas. 
}
  \item{working.cor}{
Assume working correlation structure. Could be ``unstructured" (not recommended) ``independence", ``exchangeable", or ``ar1".
}
  \item{result.file.name}{
Name of file to print analysis results to. 
}
}
\value{
A data frame where each row provides a vector with probes in the cluster,  the estimated effect size of exposure on the cluster, 
estimated (sandwich) standard error, p-value, the number of sites in the cluster and the chromosome in which the cluster is. 
}
\examples{
data(betas.7) ## upload methylation data
exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
rownames(covariates) <- colnames(betas.7)

### not run! to use package from Bioconductor to create a new annotation object
# source("http://bioconductor.org/biocLite.R")  
# biocLite("IlluminaHumanMethylation450k.db")
# require(IlluminaHumanMethylation450k.db)
# annot <- create.annot.triche(rownames(betas.7))

data(annot)
clusters.list <- assign.to.clusters(betas.7, annot)
GEE.results.clusters <- GEE.clusters(betas.7, clusters.list, exposure, covariates, id = colnames(betas.7), working.cor = "ex")
top.clusters.summary <- summarize.top.clusters(betas.7, covariates, exposure, id = colnames(betas.7), GEE.results.clusters, "results.tex", annot= annot)
}
