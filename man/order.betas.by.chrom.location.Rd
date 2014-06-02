\name{order.betas.by.chrom.location}
\alias{order.betas.by.chrom.location}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Organizes a matrix of methylation values to be ordered by chromosomal location
}
\description{
The matrix produces lists of ordered methylation matrices, according to location
}
\usage{
order.betas.by.chrom.location(betas, annot = NULL, annotation.file.name = NULL, return.chroms = NULL)
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
  \item{return.chroms}{
Optional list of chromosomes, if one is interested in specific chromosomes. 
}
}
\details{
Annotation: the package is adjusted to the Illumina annotation (i.e. through part of its column names) but could be easily adjusted to a an attestation data table with different column names. 

}
\value{
%%  \item{betas.by.chrom}{A list ordered by chromosome number. Each item in this list contains a matrix of methylation values for the subsets of sites from betas from the corresponding chromosome. The rows are organized according to chromosomal location. }
\item{sites.locations.by.chrom}{A list ordered by chromosome number. Each item in this list contains a matrix specifying chromosomal locations of each of the sites in the subsets of sites from betas from the corresponding chromosome. The rows are organized according to chromosomal location. }
%% ...
}
\examples{
data(betas.7) ## upload methylation data
data(annot)
dat.7.ord <- order.betas.by.chrom.location(betas.7, annot = annot)

dat.7.ord$betas.by.chrom$"7"[1:5,1:5]
dat.7.ord$sites.locations.by.chrom$"7"[1:5]
}
