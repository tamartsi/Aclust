\name{assign.to.clusters}
\alias{assign.to.clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cluster a set of probes
}
\description{
returns a list with clusters (each item is a cluster - a set of probes)
}
\usage{
assign.to.clusters(betas, annot = NULL, annotation.file.name = NULL, dist.type = "spearman", method = "single", dist.thresh = 0.5, bp.merge = 0, bp.thresh.clust = Inf, return.chroms = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
An (m by n) matrix of methylation values measured on $n$ participants in $m$ sites. 
}
  \item{annot}{
A preloaded data.table of annotation
}
  \item{annotation.file.name}{
If annot it not given, annotation.file.name can be provide the file name of the illumine annotation file and the function will load it and convert it to a data.table. If neither
annot or annotation.file.name are provided, annotation will be loaded from a Tim Triche's bioconductor R package. 
}
  \item{dist.type}{
Type of similarity distance function "spearman", "pearson" (correlation measures) or "euclid". 
}
  \item{method}{
Type of clustering function. "single", "complete" or "average". 
}
  \item{dist.thresh}{
A similarity distance threshold. Two neighboring clusters are merged to a single cluster if the similarity distance between them is above dist.thresh. 
}
  \item{bp.merge}{
A distance in chromosomal location. Any set of methylation sites within an interval smaller or equal to bp.dist will be potentially merged, depending on the similarity between sites at the ends of the interval. 
}
  \item{bp.thresh.clust}{
Optional maximum length between neighboring variables permitting to cluster them together. 
}
}
\details{
Annotation: the package is adjusted to the Illumina annotation (i.e. through part of its column names) but could be easily adjusted to a an attestation data table with different column names. 
}
\value{
A list in which each item is a vector of Illumina probe ids of methylation sites clustered together. }
\examples{
data(betas.7) ## upload methylation data
exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)

### If wants to use Tim Triche's Bioconductor annotation package
# source("http://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450k.db")
# require(IlluminaHumanMethylation450k.db)
# annot <- create.annot.triche(rownames(betas.7))

## else a version created on July 2013 is already in the package:
data(annot.7)
clusters.list <- assign.to.clusters(betas.7, annot.7)
}
