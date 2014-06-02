\name{betas.7}
\alias{betas.7}
\docType{data}
\title{
Methylation data (beta values) from chromosome 7.
}
\description{
Processed and normalized methylation data (beta values) from chromosome 7 from a study of 80 pesticides applicators. To maintain the privacy of the data,
the IDs are not related to the study IDs and go from 1 to 80, and exposure and covariates data are not provided in the package. This data set
can be used to experiment with using the clustering function and analysis pipeline, for which one can generate random exposure/covariates. 
Clustering results (but not analysis, since exposure data is not provided) obtained using this data set will match the results in the paper by Sofer et al (2013). 
}
\usage{data(betas.7)}
\format{
  The format is:
 num [1:28396, 1:80] 0.943 0.0331 0.039 0.0327 0.7211 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:28396] "cg00002531" "cg00003994" "cg00006414" "cg00006884" ...
  ..$ : chr [1:80] "subject1" "subject2" "subject3" "subject4" ...
}
\details{
Columns correspond to subjects, rows to CpG sites with Illumina ID as identifiers. 
}
\references{
Sofer et al, 2013. 
}
\examples{
data(betas.7)
betas.7[1:5,1:5]

## ordering the data by location on the chromosome (used for A-clustering):
data(annot)
dat.7.ord <- order.betas.by.chrom.location(betas.7, annot = annot)
names(dat.7.ord)
dim(dat.7.ord$betas.by.chrom$"7")
dim(dat.7.ord$sites.locations.by.chrom$"7")
dat.7.ord$sites.locations.by.chrom$"7"[1:5]
}
\keyword{datasets}
