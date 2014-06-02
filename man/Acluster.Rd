\name{Acluster}
\alias{Acluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A-clustering
}
\description{
An algorithm for clustering of adjacent clusters
}
\usage{
Acluster(ordr.vec, thresh.dist, which.clust = NULL, location.vec = NULL, max.dist = Inf, type = "single", dist.type = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ordr.vec}{
An (n by m) matrix of n measurements of m random variables (methylations). The columns are ordered. 
}
  \item{thresh.dist}{
A distance threshold. Two neighboring clusters are merged to a single cluster if the distance between them is above thresh.dist.
}
  \item{which.clust}{
A vector of initial clusters assignments of the m variables. If it is not provided, it is taken that each site is a unique cluster. 
}
  \item{location.vec}{
An m vector specifying the locations (e.g. chromosomal locations) of the variables measured in the matrix ordr.vec. 
}
  \item{max.dist}{
Optional maximum length between neighboring variables permitting to cluster them together. 
}
  \item{type}{
Type of clustering function. "single", "complete" or "average". 
}
  \item{dist.type}{
Type of distance function "spearman", "pearson", or "euclid". 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An m vector of cluster assignments of the m ordered variables. 
}
\examples{
data(betas.7)
data(annot)

dat.7.ord <- order.betas.by.chrom.location(betas.7, annot = annot)
cluster.vec <- Acluster(ordr.vec = dat.7.ord$betas.by.chrom[[1]], thresh.dist = 0.2, location.vec = dat.7.ord$sites.locations.by.chrom[[1]]$Coordinate_36, max.dist = 1000, type = "average")

cluster.vec[1:10] ## sites 1 to 6 are clustered together, but 7,8,9 are singletons. 

## alternatively, the function assign.to.clusters calls this function. 
 
}