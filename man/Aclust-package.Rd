\name{Aclust-package}
\alias{Aclust-package}
\alias{Aclust}
\docType{package}
\title{
Aclust package implements the A-clustering algorithm and pipeline for analysis of methylation data
}
\description{
An analysis of methylation data operates in a few steps. First, clustering of neighboring CpG sites to sets (some clusters will have only a single site). Then, performs GEE analysis of the
methylation regions, treated as outcomes, as affected by a single exposure. (The sandwich standard errors are used). The results are reported and optionally printed to a latex table. 
}
\details{
\tabular{ll}{
Package: \tab Aclust\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2013-07-15\cr
License: \tab What license is it under?\cr
}
First, load annotation. The proceed create a cluster list using the functions assign.to.clusters. To find association with exposure, the clusters are then analyzed using GEEs by the function
GEE.clusters. Finally, top clusters are summarized and printed to latex tables using summarize.top.clusters.
}
\author{
Tamar Sofer

Maintainer: Tamar Sofer <tsofer@hsph.havard.edu>
}
\references{
Sofer, T, Schifano, ED, Hoppin, JA, Hou^*, L and Baccarelli^* AA, 2013. ``A-clustering: A Novel Method for the Detection of Co-regulated Methylation Regions, and Regions Associated with Exposure".
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ Methylation, Regional Analysis }

\examples{
data(betas.7) ## upload methylation data
exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
rownames(covariates) <- colnames(betas.7)

data(annot)  ## load annotation created using the IlluminaHumanMethylation450k.db package on July 2013
clusters.list <- assign.to.clusters(betas.7, annot)
GEE.results.clusters <- GEE.clusters(betas.7, clusters.list, exposure, covariates, id = colnames(betas.7), working.cor = "ex")
top.clusters.summary <- summarize.top.clusters(betas.7, covariates, exposure, id = colnames(betas.7), GEE.results.clusters, "results.tex", annot= annot)
}
