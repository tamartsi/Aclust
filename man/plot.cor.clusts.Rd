\name{plot.cor.clusts}
\alias{plot.cor.clusts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots clusters on a heatmap ofa correlation matrix
}
\description{
Given a correlation matrix and a vector of cluster assignments, plots squares around clusters of at least 2 sites. 
}
\usage{
plot.cor.clusts(cor.mat, clusts, erase.non.clust = F, lwd = 2, main = "", xlab = "", labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cor.mat}{
An (m by m) correlation matrix
}
  \item{clusts}{
An $m$ vector of clusters assignments, as output by Acluster function. 
}
  \item{erase.non.clust}{
If TRUE, do not present correlations between sites that are not in the same cluster (i.e. they are just white). Default is FALSE. 
}
  \item{lwd}{
The width of the lines marking the clusters
}
  \item{main}{
Main title of the plot
}
  \item{xlab}{
Label of the X-axis
}
  \item{labels}{
Labels of the sites
}
}

\note{
A plot is shown
}


\examples{
data(betas.7)
data(annot)
dat.7.ord <- order.betas.by.chrom.location(betas.7, annot = annot)
cluster.vec <- Acluster(ordr.vec = dat.7.ord$betas.by.chrom[[1]], thresh.dist = 0.5, location.vec = dat.7.ord$sites.locations.by.chrom[[1]]$Coordinate_36, max.dist = 1000, type = "single")

c.mat <- cor(t(betas.7)[,1:15], method = "spearman")
plot.cor.clusts(cor.mat = c.mat, clusts  = cluster.vec[1:15])

}

