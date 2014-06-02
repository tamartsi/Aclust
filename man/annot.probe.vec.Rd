\name{annot.probe.vec}
\alias{annot.probe.vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate vector of Illumina probe names
}
\description{
Returns a matrix with annotations of the CpG sites represented by a given vector of probe names
}
\usage{
annot.probe.vec(probe.vec, annot = NULL, annotation.file.name = NULL, required.annotation = c("IlmnID", "Coordinate_36", "Gene_Name", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probe.vec}{
A vector of Illumina probe names
}
  \item{annot}{
A preloaded data.table of annotation
}
  \item{annotation.file.name}{
If annot it not given, annotation.file.name can be provide the file name of the illumine annotation file and the function will load it and convert it to a data.table. If neither
annot or annotation.file.name are provided, annotation will be loaded from a Tim Triche's bioconductor R package. 
}
   \item{required.annotation}{
Names of columns from Illumina annotation file, to be provided as annotation 
}
}
\value{
A matrix of requred.annotations  for each of the probes in probe.vec
}


\examples{
data(betas.7)
data(annot)
annot.sites <- annot.probe.vec(rownames(betas.7)[1:10], annot = annot)
}
