\name{summarize.top.clusters}
\alias{summarize.top.clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the top clusters
}
\description{
Returns information about the top cluster in the GEE analysis
}
\usage{
summarize.top.clusters(betas, covariates, exposure, id, clusters.GEE.results = NULL, clusters.GEE.results.file = NULL, minimum.sites = 2, top.number = 10, cutoff.fdr.pval = 0.05, cutoff.effect.size = NULL, annot = NULL, annotation.file.name = NULL, required.annotation = c("IlmnID", "Coordinate_37", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"), file.to.print.report = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
An (m by n) matrix of methylation values of $m$ methylation sites measured on $n$ individuals
}
  \item{covariates}{
Either a (n by p) vectors of adjusting covariates, ordered so that its rows corresponds to the columns of the `betas' matrix, or NULL if there are no covariates.  
}
  \item{exposure}{
A vector of size $n$ of exposure values for each individual. Same order as the covariates. 
}
  \item{id}{
A vector of size $n$ of the IDs of the individual. The IDs should match the column names of the matrix `betas'.  
}
  \item{clusters.GEE.results}{
A matrix of results from a  GEE analysis of clusters.
}
  \item{clusters.GEE.results.file}{
A file with the results from a GEE analysis of clusters (if clusters.GEE.results is not given)
}
  \item{minimum.sites}{
A minimum number of sites in a cluster to be considered in the analysis. 
}
  \item{top.number}{
The required number of top clusters to be reported. 
}
  \item{cutoff.fdr.pval}{
A significance p-value threshold (after FDR correction). Only clusters with p-value smaller than cutoff.fdr.pval will be reported. 
}
  \item{cutoff.effect.size}{
An effect size threshold. Only clusters with estimated effect size larger than cutoff.effect.size will be reported. 
}
  \item{annot}{
annotation data table. The package uses the Illumina annotation file ``illumina_450_manifest_v.1.2". 
}
  \item{annotation.file.name}{
A name of annotation file to read. By default it is not required, since one can use `annot', which is in the package. 
}
  \item{required.annotation}{
What is the annotation to be reported on each of the sites in the reported clusters? 
}
  \item{file.to.print.report}{
File name to print tables of annotation, effect sizes, and individual site analyses of sites in the top clusters. 
}
\item{print.progress}{
	Print status messages while the function progresses. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{top.clusters}{A data frame summarizing the GEE results (exposure effect estimates, etc) of the top clusters }
\item{annotation.top.clusters}{Annotations of the sites from the top clusters}
\item{individual.sites.analysis}{Results of individual sites analysis, of the sites in the top clusters.  }

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Individual sites analysis uses the robust standard error.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(betas.7)
exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
rownames(covariates) <- colnames(betas.7)

## not run! to create new annotation file from the 
## IlluminaHumanMethylation450k.db package:
## annot <- create.annot.triche(rownames(betas.7)

data(annot.7)

clusters.list <- assign.to.clusters(betas.7, annot.7)
GEE.results.clusters <- GEE.clusters(betas.7, clusters.list, exposure, covariates, id = colnames(betas.7), working.cor = "ex")
top.clusters.summary <- summarize.top.clusters(betas.7, covariates, exposure, id = colnames(betas.7), GEE.results.clusters, "results.tex", annot= annot.7)
}


