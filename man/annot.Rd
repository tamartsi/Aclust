\name{annot}
\alias{annot}
\docType{data}
\title{
Annotation for the probes on chromosome 7 in the Illumina 450K chip, created using the package IlluminaHumanMethylation450k.db from Bioconductor. 
}
\description{
a data.table with annotation for probes
}
\usage{data(annot.7)}
\format{
  A data frame with 28396 observations on the following 8 variables.
  \describe{
    \item{\code{IlmnID}}{a character vector}
    \item{\code{Infinium_Design_Type}}{a character vector}
    \item{\code{CHR}}{a character vector}
    \item{\code{Coordinate_37}}{a numeric vector}
    \item{\code{UCSC_RefGene_Name}}{a character vector}
    \item{\code{UCSC_RefGene_Group}}{a character vector}
    \item{\code{UCSC_CpG_Islands_Name}}{a character vector}
    \item{\code{Relation_to_UCSC_CpG_Island}}{a character vector}
  }
}
\examples{
data(annot.7)
head(annot.7)
}

