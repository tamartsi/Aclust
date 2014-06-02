\name{create.annot.triche}
\alias{create.annot.triche}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an annotation to a vector of probes
}
\description{
Create an annotation data table using Tim Triche's package on Bioconductor
}
\usage{
create.annot.triche(probe.vec, only.locations = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probe.vec}{
A vector of Illmina probe names 
}
  \item{only.locations}{
Are only locations of probes needed? In this case only chromosomes and SpG location will be returned. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data table with annotation information. If only.locations then chromosome and chromosomal locations are given. Otherwise, also 
infimum design type, 
gene name, functional group, associated island (or CpG resort), type of region in the resort. 
}
\example{

## not run! annotation is already in the package, created in July 2013. this is to create new annotation table. 
## data(betas.7) ## upload methylation data

## installing annotation package from Bioconductor
# source("http://bioconductor.org/biocLite.R")  
# biocLite("IlluminaHumanMethylation450k.db")
# require(IlluminaHumanMethylation450k.db)
# annot <- create.annot.triche(rownames(betas.7), only.locations = T)
}